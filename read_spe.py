#!/usr/bin/env python
"""
Read .SPE file into numpy array.

Adapted from http://wiki.scipy.org/Cookbook/Reading_SPE_files
Offsets and names taken as from SPE 3.0 File Format Specification:
ftp://ftp.princetoninstruments.com/Public/Manuals/Princeton%20Instruments/SPE%203.0%20File%20Format%20Specification.pdf
For the XML footer, see
ftp://ftp.princetoninstruments.com/public/Manuals/Princeton%20Instruments/Experiment%20XML%20Specification.pdf

Note: Use with SPE 3.0. Not backwards compatible with SPE 2.X.
"""
# TODO: make pytest modules with test_yes/no_footer.spe files
# TODO: fix docstrings to match numpy style: https://github.com/numpy/numpy/blob/master/doc/example.py
# TODO: remove pandas dependency. use builtin csv module.
# TODO: Use logging levels (and warnings.warn) instead of print.


from __future__ import absolute_import, division, print_function
import argparse
import copy
import os
import sys
import StringIO
import numpy as np
import pandas as pd


class File(object):
    """
    Handle an SPE file.
    """
    # TODO: Don't use protected keyword 'file' as class name.
    # Class-wide variables.
    _bits_per_byte = 8
    # TODO: don't hardcode number of metadata, get from user or footer if it exists
    # Assuming metadata datatype is 64-bit signed integer
    # from XML footer metadata using previous experiments with LightField.
    # Assuming metadata: time_stamp_exposure_started, time_stamp_exposure_ended, frame_tracking_number
    _num_metadata = 3
    _metadata_ntype = np.int64
    _spe_30_required_offsets = [6, 18, 34, 42, 108, 656, 658, 664, 678, 1446, 1992, 2996, 4098]
    _ntype_to_bits = {np.int8: 8, np.uint8: 8,
                      np.int16: 16, np.uint16: 16,
                      np.int32: 32, np.uint32: 32,
                      np.int64: 64, np.uint64: 64,
                      np.float32: 32, np.float64: 64}
    # Datatypes 6, 2, 1, 5 are for only SPE 2.X, not SPE 3.0.
    _datatype_to_ntype = {6: np.uint8, 3: np.uint16,
                          2: np.int16, 8: np.uint32,
                          1: np.int32, 0: np.float32,
                          5: np.float64}
    _binary_to_ntype = {"8s": np.int8, "8u": np.uint8,
                        "16s": np.int16, "16u": np.uint16,
                        "32s": np.int32, "32u": np.uint32,
                        "64s": np.int64, "64u": np.uint64,
                        "32f": np.float32, "64f": np.float64}

    def __init__(self, fname):
        """
        Initialize file.
        Open, load header and footer metadata, set current frame index.
        """
        # For online analysis, read metadata from binary header.
        # For final reductions, read more complete metadata from XML footer.
        self._fname = fname
        self._check_spe()
        self._fid = open(fname, 'rb')
        self._load_header_metadata()
        self._load_footer_metadata()
        self.current_frame_idx = 0
        return None

    # TODO: make __del__ method to close file automatically.

    def _check_spe(self):
        """
        Check that the file exists and is .spe.
        """
        if not os.path.isfile(self._fname):
            raise IOError(("File does not exist: {fname}").format(fname=self._fname))
        (fbase, fext) = os.path.splitext(self._fname)
        if fext != '.spe':
            raise IOError(("File extension not '.spe': {fname}").format(fname=self._fname))
        return None
    
    def _read_at(self, offset, size, ntype):
        """
        Seek to offset byte position then read size number of bytes in ntype format from file.
        """
        self._fid.seek(offset)
        result = np.fromfile(self._fid, ntype, int(size))
        return result

    def _load_header_metadata(self):
        """
        Load SPE metadata from binary header into a pandas dataframe
        and save as an object attribute.
        Use metadata from header for online analysis
        since XML footer does not yet exist while taking data.
        Only the fields required for SPE 3.0 files are loaded. All other fields are numpy NaN.
        See SPE 3.0 File Format Specification:
        ftp://ftp.princetoninstruments.com/Public/Manuals/Princeton%20Instruments/
        SPE%203.0%20File%20Format%20Specification.pdf
        """
        # TODO: don't remake temp file every time header metadata is needed.
        # file_header_ver and xml_footer_offset are
        # the only required header fields for SPE 3.0.
        # Header information from SPE 3.0 File Specification, Appendix A.
        # Remove comments from header CSV.
        ffmt = os.path.join(os.path.dirname(__file__), 'spe_30_header_format.csv')
        ffmt_base, ext = os.path.splitext(ffmt)
        if not os.path.isfile(ffmt):
            raise IOError("SPE 3.0 header format file does not exist: {fname}".format(fname=ffmt))
        if ext != '.csv':
            raise TypeError("SPE 3.0 header format file is not .csv: {fname}".format(fname=ffmt))
        with open(ffmt, 'rb') as fcmts:
            header_cmts = fcmts.readlines()
        header_nocmts = []
        for line in header_cmts:
            if line.startswith('#'):
                continue
            else:
                header_nocmts.append(line)
        header_nocmts_str = StringIO.StringIO(''.join(header_nocmts))
        self.header_metadata = pd.read_csv(header_nocmts_str, sep=',')
        # TODO: Efficiently read values and create column following
        # http://pandas.pydata.org/pandas-docs/version/0.13.1/cookbook.html
        # Index values by offset byte position.
        offset_to_value = {}
        for idx in xrange(len(self.header_metadata)):
            offset = int(self.header_metadata["Offset"][idx])
            try:
                size = (self.header_metadata["Offset"][idx+1]
                        - self.header_metadata["Offset"][idx]
                        - 1)
            # Key error if at last value in the header
            except KeyError:
                size = 1
            ntype = File._binary_to_ntype[self.header_metadata["Binary"][idx]]
            offset_to_value[offset] = self._read_at(offset, size, ntype)
        # Store only the values for the byte offsets required of SPE 3.0 files.
        # Read only first element of these values since for files written by LightField,
        # other elements and values from offets are 0.
        nan_array = np.empty(len(self.header_metadata))
        nan_array[:] = np.nan
        self.header_metadata["Value"] = pd.DataFrame(nan_array)
        for offset in File._spe_30_required_offsets:
            tf_mask = (self.header_metadata["Offset"] == offset)
            self.header_metadata["Value"].loc[tf_mask] = offset_to_value[offset][0]
        # Check for SPE 3.0
        tf_mask = (self.header_metadata["Type_Name"] == "file_header_ver")
        version = self.header_metadata[tf_mask]["Value"].values[0]
        if version != 3:
            print(("WARNING: File is not SPE version 3.\n"
                   +" SPE version: {ver}").format(ver=version), file=sys.stderr)
        return None
    
    def _load_footer_metadata(self):
        """
        Load SPE metadata from XML footer as a string
        and save as an object attribute.
        Use metadata from footer for final reductions
        since XML footer is more complete.
        """
        tf_mask = (self.header_metadata["Type_Name"] == "XMLOffset")
        xml_offset = int(self.header_metadata[tf_mask]["Value"].values[0])
        if xml_offset == 0:
            print(("INFO: XML footer metadata is empty for:\n"
                  +" {fname}").format(fname=self._fname))
        else:
            # All XML footer metadata is contained within one line.
            # Strip anything before '<SpeFormat' or after 'SpeFormat>'
            # The byte offset for the start of the XML file can be off if the file
            # is not begun/ended correctly. Search for the beginning of the XML
            # footer 1KB before the byte offset and trim the excess.
            self._fid.seek(xml_offset - 1024)
            xml_orig = self._fid.read()
            xml_trim = copy.copy(xml_orig)
            pieces = xml_trim.partition('<SpeFormat')
            xml_trim = ''.join(pieces[1:])
            pieces = xml_trim.rpartition('SpeFormat>')
            xml_trim = ''.join(pieces[:-1])
            if xml_trim == '':
                print(("WARNING: XML footer was not partitioned correctly\n" +
                       "and may need to be reformatted."), file=sys.stderr)
                xml = xml_orig
            else:
                xml = xml_trim
            self.footer_metadata = xml
        return None


    def _get_start_offset(self):
        """
        Return offset byte position of start of all data.
        """
        # TODO: use footer metadata if it exists.
        tf_mask = (self.header_metadata["Type_Name"] == "lastvalue")
        start_offset = int(self.header_metadata[tf_mask]["Offset"].values[0] + 2)
        return start_offset

    def _get_eof_offset(self):
        """
        Return end-of-file byte position.
        """
        # TODO: use footer metadata if it exists.
        self._fid.seek(0, 2)
        eof_offset = int(self._fid.tell())
        return eof_offset

    def _get_xdim(self):
        """
        Return number of pixels along frame x-axis.
        """
        # TODO: use footer metadata if it exists
        tf_mask = (self.header_metadata["Type_Name"] == "xdim")
        xdim = int(self.header_metadata[tf_mask]["Value"].values[0])
        return xdim

    def _get_ydim(self):
        """
        Return number of pixels along frame y-axis.
        """
        # TODO: use footer metadata if it exists
        tf_mask = (self.header_metadata["Type_Name"] == "ydim")
        ydim = int(self.header_metadata[tf_mask]["Value"].values[0])
        return ydim

    def _get_pixels_per_frame(self):
        """
        Return number of pixels per frame.
        """
        # TODO: use footer metadata if it exists.
        xdim = self._get_xdim()
        ydim = self._get_ydim()
        pixels_per_frame = int(xdim * ydim)
        return pixels_per_frame

    def _get_pixel_ntype(self):
        """
        Return pixel binary data type as numpy type.
        """
        # TODO: use footer metadata if it exists.
        tf_mask = (self.header_metadata["Type_Name"] == "datatype")
        pixel_datatype = self.header_metadata[tf_mask]["Value"].values[0]
        pixel_ntype = File._datatype_to_ntype[pixel_datatype]
        return pixel_ntype

    def _get_bytes_per_frame(self):
        """
        Return number of bytes per frame.
        """
        # TODO: use footer metadata if it exists.
        # Infer frame size.
        # From SPE 3.0 File Format Specification, Ch 1 (with clarifications):
        # bytes_per_frame = pixels_per_frame * bits_per_pixel / (8 bits per byte)
        pixels_per_frame = self._get_pixels_per_frame()
        pixel_ntype = self._get_pixel_ntype()
        bits_per_pixel = File._ntype_to_bits[pixel_ntype]
        bytes_per_frame = int(pixels_per_frame * (bits_per_pixel / File._bits_per_byte))
        return bytes_per_frame

    def _get_bytes_per_metadata_elt(self):
        """
        Return number of bytes per element of metadata.
        """
        # TODO: use footer metadata if it exists.
        # From SPE 3.0 File Format Specification, Ch 1 (with clarifications):
        # bytes_per_metadata_elt = 8 bytes per metadata element
        #   metadata element includes time stamps, frame tracking number, etc with 8 bytes each.
        bits_per_metadata_elt = File._ntype_to_bits[File._metadata_ntype]
        bytes_per_metadata_elt = int(bits_per_metadata_elt / File._bits_per_byte)
        return bytes_per_metadata_elt

    def _get_bytes_per_metadata_set(self):
        """
        Return number of bytes per set of metadata elements.
        """
        bytes_per_metadata_elt = self._get_bytes_per_metadata_elt()
        bytes_per_metadata_set = int(File._num_metadata * bytes_per_metadata_elt)
        return bytes_per_metadata_set

    def _get_bytes_per_stride(self):
        """
        Return number of bytes per frame + per-frame metadata.
        Equivalent to the number of bytes to move to the beginning of the next frame.
        """
        bytes_per_frame = self._get_bytes_per_frame()
        bytes_per_metadata_set = self._get_bytes_per_metadata_set()
        bytes_per_stride = int(bytes_per_frame + bytes_per_metadata_set)
        return bytes_per_stride
        
    def get_num_frames(self):
        """
        Return number of frames currently in an SPE file.
        """
        # TODO: use footer metadata if it exists.
        # Infer the number of frames that have been taken using the file size in bytes.
        # NumFrames from the binary header metadata is the 
        # number of frames typed into LightField that will potentially be taken,
        # not the number of frames that have already been taken and are in the file being read.
        # In case the file is currently being written to by LightField
        # when the file is being read by Python, count only an integer number of frames.
        # Allow negative indexes using mod.
        start_offset = self._get_start_offset()
        bytes_per_stride = self._get_bytes_per_stride()
        eof_offset = self._get_eof_offset()
        num_frames = int((eof_offset - start_offset) // bytes_per_stride)
        return num_frames
                
    def get_frame(self, frame_idx):
        """
        Return a frame and per-frame metadata from the file.
        Frame is returned as a numpy 2D array.
        Time stamp metadata is returned as Python datetime object.
        frame_idx argument is python indexed: 0 is first frame.
        """
        # See SPE 3.0 File Format Specification:
        # ftp://ftp.princetoninstruments.com/Public/Manuals/Princeton%20Instruments/
        # SPE%203.0%20File%20Format%20Specification.pdf
        # If XML footer metadata exists (i.e. for final reductions).
        if hasattr(self, 'footer_metadata'):
            # TODO: complete as below
            pass
        # Else use binary header metadata (i.e. for online analysis).
        # else:
        # Get the number of frames currently in the file.
        # Update the index position of the frame last read.
        num_frames = self.get_num_frames()
        self.current_frame_idx = int(frame_idx % num_frames)
        # Infer frame and per-frame metadata offsets.
        start_offset = self._get_start_offset()
        bytes_per_stride = self._get_bytes_per_stride()
        frame_offset = start_offset + (self.current_frame_idx * bytes_per_stride)
        bytes_per_frame = self._get_bytes_per_frame()
        metadata_offset = frame_offset + bytes_per_frame
        bytes_per_metadata_elt = self._get_bytes_per_metadata_elt()
        # Read frame, metadata. Format metadata timestamps to be absolute time, UTC.
        # Time_stamps from the ProEM's internal timer-counter card are in 1E6 ticks per second.
        # Ticks per second from XML footer metadata using previous LightField experiments:
        # 1 tick = 1 microsecond ; 1E6 ticks per second.
        # 0 ticks is when "Acquire" was first clicked on LightField.
        pixels_per_frame = self._get_pixels_per_frame()
        pixel_ntype = self._get_pixel_ntype()
        frame = self._read_at(frame_offset, pixels_per_frame, pixel_ntype)
        xdim = self._get_xdim()
        ydim = self._get_ydim()
        frame = frame.reshape((ydim, xdim))
        mtsexpstart_offset = metadata_offset
        mtsexpend_offset = mtsexpstart_offset + bytes_per_metadata_elt
        mftracknum_offset = mtsexpend_offset + bytes_per_metadata_elt
        metadata = {}
        mtsexpstart = self._read_at(mtsexpstart_offset, 1, File._metadata_ntype)[0]
        mtsexpend   = self._read_at(mtsexpend_offset, 1, File._metadata_ntype)[0]
        mftracknum  = self._read_at(mftracknum_offset, 1, File._metadata_ntype)[0]
        metadata["time_stamp_exposure_started"] = mtsexpstart
        metadata["time_stamp_exposure_ended"] = mtsexpend
        metadata["frame_tracking_number"] = mftracknum
        return (frame, metadata)

    # # TODO: make generator. for now use get_frame
    # def get_frames(self, frame_idx_list):
    #     """
    #     Yield a frame and per-frame metadata from the file.
    #     Return a frame and per-frame metadata from the file.
    #     Frame is returned as a numpy 2D array.
    #     Time stamp metadata is returned as Python datetime object.
    #     frame_list argument is python indexed: 0 is first frame.
    #     """
    #     # get_num_frames()
    #     # self.current_frame_idx
    #     for fnum in frame_idx_list:
    #         print(fnum)
    #     return None

    def close(self):
        """
        Close file.
        """
        self._fid.close()
        return None

def main(args):
    """
    Read a numbered frame from the SPE file.
    Show a plot and print the metadata.
    """
    fid = File(args.fname)
    (frame, metadata) = fid.get_frame(args.frame_idx)
    fid.close()
    return (frame, metadata)
            
if __name__ == "__main__":
    # TODO: have defaults for metadata
    fname_default = "test_yes_footer.spe"
    frame_idx_default = -1
    parser = argparse.ArgumentParser(description="Read a SPE file and return ndarray frame and dict metadata variables.")
    parser.add_argument("--fname",
                        default=fname_default,
                        help=("Path to SPE file. "
                              +"Default: {default}".format(default=fname_default)))
    parser.add_argument("--frame_idx",
                        default=frame_idx_default,
                        help=("Frame index to read in. First frame is 0. Last frame is -1. "
                              +"Default: {default}".format(default=frame_idx_default)))
    parser.add_argument("--verbose",
                        "-v",
                        action='store_true',
                        help=("Print 'INFO:' messages to stdout."))
    args = parser.parse_args()
    if args.verbose:
        print("INFO: Arguments:")
        for arg in args.__dict__:
            print(arg, args.__dict__[arg])
    (frame, metadata) = main(args)
    
