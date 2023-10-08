# SpectrumAnalysis is a program of automatic analysis of bone FTIRM spectrum.
# Copyright (C) 2017  GARDEGARONT Marc
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from array import array
import numpy as np
import matplotlib.pyplot as plt

def spload2(filename):
    DSet2DC1DIBlock               =  120;
    HistoryRecordBlock            =  121;
    InstrHdrHistoryRecordBlock    =  122;
    InstrumentHeaderBlock         =  123;
    IRInstrumentHeaderBlock       =  124;
    UVInstrumentHeaderBlock       =  125;
    FLInstrumentHeaderBlock       =  126;
    # Data member IDs
    DataSetDataTypeMember              =  -29839;
    DataSetAbscissaRangeMember         =  -29838;
    DataSetOrdinateRangeMember         =  -29837;
    DataSetIntervalMember              =  -29836;
    DataSetNumPointsMember             =  -29835;
    DataSetSamplingMethodMember        =  -29834;
    DataSetXAxisLabelMember            =  -29833;
    DataSetYAxisLabelMember            =  -29832;
    DataSetXAxisUnitTypeMember         =  -29831;
    DataSetYAxisUnitTypeMember         =  -29830;
    DataSetFileTypeMember              =  -29829;
    DataSetDataMember                  =  -29828;
    DataSetNameMember                  =  -29827;
    DataSetChecksumMember              =  -29826;
    DataSetHistoryRecordMember         =  -29825;
    DataSetInvalidRegionMember         =  -29824;
    DataSetAliasMember                 =  -29823;
    DataSetVXIRAccyHdrMember           =  -29822;
    DataSetVXIRQualHdrMember           =  -29821;
    DataSetEventMarkersMember          =  -29820;
    # Type code IDs
    ShortType               = 29999;
    UShortType              = 29998;
    IntType                 = 29997;
    UIntType                = 29996;
    LongType                = 29995;
    BoolType                = 29988;
    CharType                = 29987;
    CvCoOrdPointType        = 29986;
    StdFontType             = 29985;
    CvCoOrdDimensionType    = 29984;
    CvCoOrdRectangleType    = 29983;
    RGBColorType            = 29982;
    CvCoOrdRangeType        = 29981;
    DoubleType              = 29980;
    CvCoOrdType             = 29979;
    ULongType               = 29978;
    PeakType                = 29977;
    CoOrdType               = 29976;
    RangeType               = 29975;
    CvCoOrdArrayType        = 29974;
    EnumType                = 29973;
    LogFontType             = 29972;


    fid = open(filename,'rb')
    if fid == -1:
        error('Cannot open the file.')
        return

    # a = array('u')
    # Fixed file header of signature and description
    signature = ''.join([chr(c) for c in np.fromfile(fid, np.uint8, count=4)])  # ORIGINAL : char(fread(fid, 4, 'uchar')')

    if signature != 'PEPE':

        error('This is not a PerkinElmer block structured file.')
        return
    description = ''.join([chr(c) for c in np.fromfile(fid, np.uint8, count=40)])  # ORIGINAL : char(fread(fid, 40, 'uchar')')

    # Initialize a variable so we can tell if we have read it.
    xLen = np.int32(0)

    # The rest of the file is a list of blocks
    while True:
        try:
            blockID = np.fromfile(fid, np.int16, count=1)[0]
            blockSize = np.fromfile(fid, np.int32, count=1)[0]
        except:
            break

        # feof does not go true until after the read has failed.
        #if feof(fid):
        #    break

        if blockID == DSet2DC1DIBlock:
            pass
            # Wrapper block.  Read nothing.

        elif blockID == DataSetAbscissaRangeMember:
            innerCode = np.fromfile(fid, np.int16, count=1)[0]
            #_ASSERTE(CvCoOrdRangeType == nInnerCode);
            x0 = np.fromfile(fid, np.float, count=1)[0]
            xEnd = np.fromfile(fid, np.float, count=1)[0]

        elif blockID == DataSetIntervalMember:
            innerCode = np.fromfile(fid, np.int16, count=1)[0]
            xDelta = np.fromfile(fid, np.float, count=1)[0]

        elif blockID == DataSetNumPointsMember:
            innerCode = np.fromfile(fid, np.int16, count=1)[0]
            xLen = np.fromfile(fid, np.int32, count=1)[0]

        elif blockID == DataSetXAxisLabelMember:
            innerCode = np.fromfile(fid, np.int16, count=1)[0]
            len = np.fromfile(fid, np.int16, count=1)[0]
            xLabel = ''.join([chr(c) for c in np.fromfile(fid, np.uint8, count=len)])

        elif blockID == DataSetYAxisLabelMember:
            innerCode = np.fromfile(fid, np.int16, count=1)[0]
            len = np.fromfile(fid, np.int16, count=1)[0]
            yLabel = ''.join([chr(c) for c in np.fromfile(fid, np.uint8, count=len)])

        elif blockID == DataSetAliasMember:
            innerCode = np.fromfile(fid, np.int16, count=1)[0]
            len = np.fromfile(fid, np.int16, count=1)[0]
            alias = ''.join([chr(c) for c in np.fromfile(fid, np.uint8, count=len)])

        elif blockID == DataSetNameMember:
            innerCode = np.fromfile(fid, np.int16, count=1)[0]
            len = np.fromfile(fid, np.int16, count=1)[0]
            originalName = ''.join([chr(c) for c in np.fromfile(fid, np.uint8, count=len)])

        elif blockID == DataSetDataMember:
            innerCode = np.fromfile(fid, np.int16, count=1)[0]
            len = np.fromfile(fid, np.int32, count=1)[0]
            # innerCode should be CvCoOrdArrayType
            # len should be xLen * 8
            if xLen == 0:
                xLen = len / 8
            data = np.fromfile(fid, np.float, count=xLen)

        else:
            fid.seek(blockSize, 1)
    fid.close()

    if xLen == 0:
        error('The file does not contain spectral data.')
        return

    # Expand the axes specifications into vectors
    xAxis = np.arange(x0, xEnd+xDelta, xDelta)

    return [data, xAxis]
"""
    # Return the other details as name,value pairs
    misc(1,:) = {'xLabel', xLabel};
    misc(2,:) = {'yLabel', yLabel};
    misc(3,:) = {'alias', alias};
    misc(4,:) = {'original name', originalName};
"""
