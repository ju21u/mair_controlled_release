from netCDF4 import Dataset as Dataset
import numpy.ma as ma
import glob
import shutil
import argparse
import os


def create_wrfchemi(input_dir: str, output_dir: str, new_date: str, new_time: str = "00:00:00"):
    print(
        f"input_dir: {input_dir}, output_dir: {output_dir}, new_date: {new_date}, new_time: {new_time}"
    )
    list_of_wrfchemi_files = glob.glob(f"{input_dir}/wrfchemi*")

    for current_wrfchemi in list_of_wrfchemi_files:
        # filename format: wrfchemi_d01_2017-09-06_12_00_00
        filename = os.path.basename(current_wrfchemi)
        # extract non date portion of filename
        prefix = filename[0:13]
        suffix = "_" + new_time  # if you want to specify the time
        new_filename = prefix + new_date + suffix
        output_wrfchemi = f"{output_dir}/{new_filename}"
        shutil.copyfile(src=current_wrfchemi, dst=output_wrfchemi)

        ds = Dataset(output_wrfchemi, "r+")
        tstep = 0

        new_datetime = new_date + "_" + new_time
        new_array = [bytes(i, "ascii") for i in new_datetime]

        ds["Times"][tstep] = new_array

        ds.close()


def main():
    parser = argparse.ArgumentParser(description="Create WRFChemi files")

    parser.add_argument("--input-dir", required=True, help="directory with initial wrfchemi")
    parser.add_argument("--output-dir", required=True, help="directory to copy modified wrfchemi")
    parser.add_argument(
        "--new-date", required=True, help="new date for wrfchemi in form YYYY-MM-DD"
    )
    parser.add_argument("--new-time", help="new date for wrfchemi in form HH:MM:SS")
    args = parser.parse_args()

    create_wrfchemi(args.input_dir, args.output_dir, args.new_date, args.new_time)


if __name__ == "__main__":
    main()
