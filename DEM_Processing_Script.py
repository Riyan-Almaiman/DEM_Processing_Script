import os
import time
import datetime
import concurrent.futures
from threading import Thread
import arcpy
import subprocess
from osgeo import gdal
import rasterio
import numpy as np
from shapely.geometry import Point
import geopandas as gpd
import multiprocessing
# Constants
MAX_MOSAIC_THREADS = 10
MAX_MINUS_THREADS = 5
osgeo4w_path = r"C:\OSGeo4W\OSGeo4W.bat"  
thresholds = [10, 40]

def calculate_tri(input_dem, output_tri):
    dem = gdal.Open(input_dem)
    if not dem:
        raise IOError("Could not open DEM file")
    gdal.DEMProcessing(output_tri, dem, 'TRI')
    dem = None
    

def create_mosaic_arcpy(tifs, output_folder, raster_dataset_name, pixel_type='32_BIT_FLOAT', bands=1):
    if not os.path.exists(os.path.join(output_folder, raster_dataset_name)):
        print(f"Starting Mosaic for folder: {raster_dataset_name}")
        start_time = time.time()
        arcpy.env.workspace = output_folder
        arcpy.MosaicToNewRaster_management(
            input_rasters=tifs,
            output_location=output_folder,
            raster_dataset_name_with_extension=raster_dataset_name,
            pixel_type=pixel_type,
            number_of_bands=bands
        )
        end_time = time.time()
        print(f"Mosaic creation for {raster_dataset_name} completed in {(end_time - start_time) / 60:.2f} minutes.")

def process_folder_arcpy(folder, base_folder, destination_selection, tri_directory, shp_base_directory):
    DEM_types = ['DTM', 'DSM']
    folder_path = os.path.join(base_folder, folder)
    if os.path.isdir(folder_path):
        for DEM_type in DEM_types:
            tif_files = [os.path.join(folder_path, "data", DEM_type.lower(), "data", f)
                         for f in os.listdir(os.path.join(folder_path, "data", DEM_type.lower(), "data"))
                         if f.endswith('.tif')]
            if tif_files:
                output_folder = os.path.join(destination_selection, f"DEM_{folder.upper()}", f"DEM_{folder.upper()}_{DEM_type}")
                os.makedirs(output_folder, exist_ok=True) 
                raster_dataset_name = f"{folder.upper()}_{DEM_type}.tif"
                create_mosaic_arcpy(tif_files, output_folder, raster_dataset_name)
                if DEM_type == 'DTM':
                    dtm_folder_path = os.path.join(output_folder, raster_dataset_name)
                    createTRI(dtm_folder_path, tri_directory, shp_base_directory)

def process_folders_multiprocessing(base_folder, destination_selection, tri_directory, shp_base_directory):
    folders = [folder for folder in os.listdir(base_folder) if os.path.isdir(os.path.join(base_folder, folder))]
    with multiprocessing.Pool(processes=MAX_MOSAIC_THREADS) as pool:
        pool.starmap(process_folder_arcpy, [(folder, base_folder, destination_selection, tri_directory, shp_base_directory) for folder in folders])


def find_dsm_dtm_folders(geocell_folder):
    dsm_folder = os.path.join(geocell_folder, "data/dsm/data")
    dtm_folder = os.path.join(geocell_folder, "data/dtm/data")
    return dsm_folder, dtm_folder

def run_gdal_command(dsm_folder, dtm_folder, output_folder_name, destination_directory):
    output_dir = os.path.join(destination_directory, output_folder_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        start_time = time.time()
        for dsm_file in os.listdir(dsm_folder):
            dtm_file = dsm_file.replace('dsm', 'dtm')
            if os.path.exists(os.path.join(dtm_folder, dtm_file)):
                dsm_path = os.path.join(dsm_folder, dsm_file)
                dtm_path = os.path.join(dtm_folder, dtm_file)
                output_path = os.path.join(output_dir, dsm_file)
                gdal_command = f"gdal_calc.py -A {dsm_path} -B {dtm_path} --outfile {output_path} --calc=\"A-B\""
                with subprocess.Popen(osgeo4w_path, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True, shell=True) as proc:
                    proc.stdin.write(gdal_command + '\n')
                    proc.stdin.close()
                    proc.wait()
        end_time = time.time()
        print(f"Minus operation for {output_folder_name} completed in {(end_time - start_time) / 60:.2f} minutes.")

def create_minus_parallel(raw_folders, destination_directory):
    folder_names = [d for d in os.listdir(raw_folders) if os.path.isdir(os.path.join(raw_folders, d))]
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_MINUS_THREADS) as executor:
        futures = [executor.submit(run_gdal_command, find_dsm_dtm_folders(os.path.join(raw_folders, folder))[0], find_dsm_dtm_folders(os.path.join(raw_folders, folder))[1], folder, destination_directory) for folder in folder_names]
        concurrent.futures.wait(futures)

def createTRI(mosaic_path, tri_directory, shp_base_directory):
    dem_folders = [d for d in os.listdir(mosaic_path) if os.path.isdir(os.path.join(mosaic_path, d)) and d.startswith('DEM_')]
    for dem_folder in dem_folders:
        dtm_folder_path = os.path.join(mosaic_path, dem_folder, f"{dem_folder}_DTM")
        if os.path.exists(dtm_folder_path):
            tif_files = [f for f in os.listdir(dtm_folder_path) if f.endswith('.tif')]
            for dem_file in tif_files:
                dem_path = os.path.join(dtm_folder_path, dem_file)
                dem_base_name = os.path.splitext(dem_file)[0]
                tri_output_path = os.path.join(tri_directory, f"{dem_base_name}_TRI.tif")
                if not os.path.exists(tri_output_path):
                    try:
                        calculate_tri(dem_path, tri_output_path)
                        print(f"Created TRI: {tri_output_path}")
                    except Exception as e:
                        print(f"Error generating TRI: {e}")
                shp_directories = {threshold: os.path.join(shp_base_directory, f"threshold_{threshold}") for threshold in thresholds}
                for threshold, shp_directory in shp_directories.items():
                    os.makedirs(shp_directory, exist_ok=True)
                    shp_output_path = os.path.join(shp_directory, f"{dem_base_name}_wells_{threshold}.shp")
                    if not os.path.exists(shp_output_path):
                        with rasterio.open(tri_output_path) as src:
                            tri_data = src.read(1)
                            high_ruggedness_pixels = tri_data > threshold
                            rows, cols = np.where(high_ruggedness_pixels)
                            if rows.size and cols.size:
                                xs, ys = rasterio.transform.xy(src.transform, rows, cols, offset='center')
                                points = [Point(x, y) for x, y in zip(xs, ys)]
                                gdf = gpd.GeoDataFrame(geometry=points, crs=src.crs)
                                gdf.to_file(shp_output_path)

def create_shp_files(tri_output_path, shp_base_directory, dem_base_name):
    shp_directories = {threshold: os.path.join(shp_base_directory, f"threshold_{threshold}") for threshold in thresholds}
    for threshold, shp_directory in shp_directories.items():
        os.makedirs(shp_directory, exist_ok=True)
        shp_output_path = os.path.join(shp_directory, f"{dem_base_name}_wells_{threshold}.shp")
        if not os.path.exists(shp_output_path):
            with rasterio.open(tri_output_path) as src:
                tri_data = src.read(1)
                high_ruggedness_pixels = tri_data > threshold
                rows, cols = np.where(high_ruggedness_pixels)
                if rows.size and cols.size:
                    xs, ys = rasterio.transform.xy(src.transform, rows, cols, offset='center')
                    points = [Point(x, y) for x, y in zip(xs, ys)]
                    gdf = gpd.GeoDataFrame(geometry=points, crs=src.crs)
                    gdf.to_file(shp_output_path)


def check_new_folders_and_process(base_dir, check_interval=1800):
    processed_folders = set()
    while True:
        current_time = datetime.datetime.now()
        print(f"Checking for new folders at {current_time.strftime('%Y-%m-%d %H:%M:%S')}")

        for folder in os.listdir(base_dir):
            folder_path = os.path.join(base_dir, folder)
            raw_folder_path = os.path.join(folder_path, 'raw')
            done_file_path = os.path.join(folder_path, 'done.txt')

            if os.path.isdir(raw_folder_path) and os.path.isfile(done_file_path) and folder not in processed_folders:
                print(f"Processing new folder: {folder}")
                processed_folders.add(folder)
                mosaic_path, minus_path, TRI_Path, TRI_SHP_Path = setup_directories(folder_path)

                mosaic_thread = Thread(target=process_folders_multiprocessing, args=(raw_folder_path, mosaic_path, TRI_Path, TRI_SHP_Path))
                minus_thread = Thread(target=create_minus_parallel, args=(raw_folder_path, minus_path))
                mosaic_thread.start()
                minus_thread.start()

                mosaic_thread.join()
                minus_thread.join()
                createTRI(mosaic_path, TRI_Path, TRI_SHP_Path)


                os.rename(done_file_path, os.path.join(folder_path, 'processing_done.txt'))
                print(f"Finished processing folder: {folder}, renamed done file to processing_done.txt")

        time.sleep(check_interval)



def setup_directories(folder_path):
    mosaic_path = os.path.join(folder_path, 'Mosaic')
    minus_path = os.path.join(folder_path, 'Minus')
    TRI_Path = os.path.join(folder_path, 'TRI')
    TRI_SHP_Path = os.path.join(folder_path, 'TRI_SHP')
    for path in [mosaic_path, minus_path, TRI_Path, TRI_SHP_Path]:
        if not os.path.exists(path):
            os.mkdir(path)
    return mosaic_path, minus_path, TRI_Path, TRI_SHP_Path

if __name__ == "__main__":
    base_directory = r""
    check_new_folders_and_process(base_directory)
