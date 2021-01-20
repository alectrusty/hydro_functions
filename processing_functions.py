"""
# processing_functions.py
"""
import os
import whitebox
import richdem as rd
from osgeo import gdal
import geopandas as gpd
import math
import numpy
from tqdm import tqdm


def clear():
    os.system('cls')

def create_output_structure(output_path):
    """
    Generate output folder structure.
    """
    try:
        # create an output folder for the OpenHydro results
        # oh_output_path = os.path.join(output_path, 'OpenHydro_Output')
        oh_output_path = output_path

        # populate sub-folders for tool output
        os.makedirs(os.path.join(oh_output_path, '00_Huc10'))
        os.makedirs(os.path.join(oh_output_path, '01_Geomorphons'))
        os.makedirs(os.path.join(oh_output_path, '02_Topographic_Openness'))
        os.makedirs(os.path.join(oh_output_path, '03_Planform_Curvature'))
        os.makedirs(os.path.join(oh_output_path, '04_Smoothed'))
        os.makedirs(os.path.join(oh_output_path, '05_Fill'))
        os.makedirs(os.path.join(oh_output_path, '06_Sinks'))
        os.makedirs(os.path.join(oh_output_path, '07_Flow_Direction'))
        os.makedirs(os.path.join(oh_output_path, '08_Flow_Accumulation'))
        os.makedirs(os.path.join(oh_output_path, '09_Stream_Definition'))
        os.makedirs(os.path.join(oh_output_path, '10_Stream_Polylines'))

    except OSError as e:
        print("Output folder already exists, moving on.")
        pass

def prepare_dem(input_huc_file,
                input_dem_file,
                output_path,
                huc_name,
                epsg,
                buffer_size):
    """
    Adds a user defined spatial buffer to a given HUC10 shapefile,
    then clips the DTM to boundary of the provided polygon.
    """

    # set up OH output path
    # oh_output_path = os.path.join(output_path, 'OpenHydro_Output')
    oh_output_path = output_path

    # TODO: Add a check for the units of the dataset -- make sure it's in meters
    # read the HUC10 feature dataset
    huc_shape = gpd.read_file(input_huc_file).to_crs(epsg)
    
    # buffer by the provided distance
    huc_buffed = huc_shape.buffer(buffer_size)
    
    # write out the buffered huc shapefile to the output path
    clipped_huc_file = os.path.join(oh_output_path, '00_Huc10/' + huc_name + "_" + str(buffer_size) + '_buffer.shp')
    huc_buffed.to_file(clipped_huc_file)
    
    # clip the dem to that buffered shapefile using a gdal warp call
    clipped_dem_file = os.path.join(oh_output_path, '00_Huc10/' + huc_name + "_" + str(buffer_size) + '_dem.tif')
    gdal.Warp(clipped_dem_file, input_dem_file, cutlineDSName=clipped_huc_file, dstSRS=epsg)
    

def flow_accumulation_workflow(input_dem_file,
                               output_path,
                               huc_name,
                               filter_size,
                               flowdir):
    """Runs a full flow accumulation workflow.
    Parameters:
    input_dem_file: A path to an input DEM
    oh_output_path: A path to the folder where all OpenHydro output goes.
    Outputs:
    smoothed_dem_file
    filled_dem_file
    sinks_dem_file
    flowdir_dem_file
    
    Returns: Printed processing updates.
   """
    

    # print("Starting flow accumulation workflow...")
    wbt = whitebox.WhiteboxTools()
    wbt.verbose = False

    # set output folder path
    # TODO: Check if the OpenHydro folder exists, if not, create one.
    # oh_output_path = os.path.join(output_path, 'OpenHydro_Output')
    oh_output_path = output_path
        
    # set output file locations
    smoothed_dem_file = os.path.join(oh_output_path, '04_Smoothed/' + huc_name + '_smoothed.tif')
    filled_dem_file = os.path.join(oh_output_path, '05_Fill/' + huc_name + '_filled.tif')
    sinks_dem_file = os.path.join(oh_output_path, '06_Sinks/' + huc_name + '_sinks.tif')
    #flowdir_dem_file = os.path.join(oh_output_path, '06_Flow_Direction/' + huc_name + '_flowdir.tif')
    
    # run feature preserving
    # print("Smoothing DEM...")
    wbt.feature_preserving_smoothing(input_dem_file, smoothed_dem_file, filter=filter_size)
    
    # run breach depressions -- i.e. 'filled' DEM
    # print("Breaching depressions...")
    wbt.breach_depressions(smoothed_dem_file, filled_dem_file)

    # run sinks
    # print("Calculating sinks...")
    wbt.sink(smoothed_dem_file, sinks_dem_file)
    
    # run flow direction
    if flowdir == 'd-inf':
        print("Starting flow accumulation using the D-Infinity method...")
        dinf_flac_dem_file = os.path.join(oh_output_path, '08_Flow_Accumulation/' + huc_name + '_d_inf_flac.tif')
        wbt.d_inf_flow_accumulation(filled_dem_file, dinf_flac_dem_file)
    elif flowdir == 'd8':
        print("Starting flow accumulation using the D-8 method...")
        d8_flac_dem_file = os.path.join(oh_output_path, '08_Flow_Accumulation/' + huc_name + '_d8_flac.tif')
        wbt.d8_flow_accumulation(filled_dem_file, d8_flac_dem_file)


def stream_extraction_workflow(filled_dem_file,
                               flac_dem_file,
                               output_path,
                               huc_name,
                               channel_threshold,
                               zero_background=False):

    """Stream Extraction Workflow"""

    # print("Starting stream extraction workflow...")
    wbt = whitebox.WhiteboxTools()
    wbt.verbose = False

    # set output folder path
    # oh_output_path = os.path.join(output_path, 'OpenHydro_Output')
    oh_output_path = output_path

    # set output file locations
    stream_def_raster = os.path.join(oh_output_path, '09_Stream_Definition/' + huc_name + "_stream_definition.tif")
    d8_pointer_raster = os.path.join(oh_output_path, '07_Flow_Direction/' + huc_name + "_d8_pointer.tif")
    stream_def_polygon = os.path.join(oh_output_path, '10_Stream_Polylines/' + huc_name + "_stream_polylines.shp")
    # stream_link_raster = os.path.join(oh_output_path, "09_Stream_Definition/" + huc_name + "_stream_link.tif")
    # stream_length_raster = os.path.join(oh_output_path, "09_Stream_Definition/" + huc_name + "_stream_length.tif")

    # extract streams given a flow accumulation raster
    # print("Extracting stream channels given threshold: " + str(channel_threshold) + " ...")
    wbt.extract_streams(flac_dem_file, stream_def_raster, channel_threshold, zero_background)

    # create a d8 pointer file (mandatory for raster_streams_to_vector)
    # print("Creating a d8 pointer file...")
    wbt.d8_pointer(filled_dem_file, d8_pointer_raster)

    # identify stream link
    # wbt.stream_link_identifier(d8_pointer_raster, stream_def_raster, stream_link_raster)
    # wbt.stream_link_length(d8_pointer_raster, stream_link_raster, stream_length_raster)

    # remove short streams from stream definition raster
    # print("Removing short streams shorter than " + str(min_seg_length) + " map units...")
    # wbt.remove_short_streams(d8_pointer_raster, stream_length_raster,  stream_def_raster, min_length=min_seg_length)

    # convert the raster stream definition file to vector
    # print("Converting stream definition raster to vector...")
    wbt.raster_streams_to_vector(stream_def_raster, d8_pointer_raster, stream_def_polygon)

    print("Stream extraction workflow complete.")


def planform_curvature(smoothed_dem_file, output_path):
    """Calculates planform curvature for a given DEM and writes it out to a specified location"""

    # load the dem
    rd_dem = rd.LoadGDAL(smoothed_dem_file)

    # calculates planform curvature
    planform_curv = rd.TerrainAttribute(rd_dem, attrib='planform_curvature')

    # saves the planform curvature raster
    rd.SaveGDAL(output_path, planform_curv)
    
    
   def retile_raster(input_tif, input_path, output_path, tile_size_x, tile_size_y):
    """Retiles a given raster to a specified dimension"""
    # get relative path
    parent_folder = os.path.dirname(os.path.realpath(__file__))
    gdal_translate_path = os.path.join(parent_folder, 'Python3/Lib/site-packages/osgeo/gdal_translate.exe')

    # subtile prefix
    output_filename = 'tile_'

    # open input tif and get relevant information
    ds = gdal.Open(input_path + input_tif)
    band = ds.GetRasterBand(1)
    xsize = band.XSize
    ysize = band.YSize

    # call the specific GDAL translate script
    # "c"
    for i in range(0, xsize, tile_size_x):
        for j in range(0, ysize, tile_size_y):
            com_string = gdal_translate_path + " -of GTIFF -srcwin " + str(
                i) + ", " + str(j) + ", " + str(
                tile_size_x) + ", " + str(tile_size_y) + " " + str(input_path) + str(input_tif) + " " + str(
                output_path) + str(output_filename) + str(i) + "_" + str(j) + ".tif"
            # print(com_string)
            os.system(com_string)


def make_processing_grid(outputGridfn,xmin,xmax,ymin,ymax,gridHeight,gridWidth):

    # convert sys.argv to float
    xmin = float(xmin)
    xmax = float(xmax)
    ymin = float(ymin)
    ymax = float(ymax)
    gridWidth = float(gridWidth)
    gridHeight = float(gridHeight)

    # get rows
    rows = ceil((ymax-ymin)/gridHeight)
    # get columns
    cols = ceil((xmax-xmin)/gridWidth)

    # start grid cell envelope
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax-gridHeight

    # create output file
    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outputGridfn):
        os.remove(outputGridfn)
    outDataSource = outDriver.CreateDataSource(outputGridfn)
    outLayer = outDataSource.CreateLayer(outputGridfn,geom_type=ogr.wkbPolygon )
    featureDefn = outLayer.GetLayerDefn()

    # create grid cells
    countcols = 0
    while countcols < cols:
        countcols += 1

        # reset envelope for rows
        ringYtop = ringYtopOrigin
        ringYbottom =ringYbottomOrigin
        countrows = 0

        while countrows < rows:
            countrows += 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            # add new geom to layer
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(poly)
            outLayer.CreateFeature(outFeature)
            outFeature.Destroy

            # new envelope for next poly
            ringYtop = ringYtop - gridHeight
            ringYbottom = ringYbottom - gridHeight

        # new envelope for next poly
        ringXleftOrigin = ringXleftOrigin + gridWidth
        ringXrightOrigin = ringXrightOrigin + gridWidth

    # Close DataSources
    outDataSource.Destroy()
