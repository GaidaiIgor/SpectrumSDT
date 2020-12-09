#!/usr/bin/env python3

import argparse
import copy
import itertools
import os
import shutil
import os.path as path
from pathlib import Path
from typing import List

from SpectrumSDTConfig import SpectrumSDTConfig
            

def parse_command_line_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Creates a folder structure for calculations with a given config")
    parser.add_argument("-c", "--config", default="spectrumsdt.config", help="Path to configuration file")
    parser.add_argument("-K", required=True, help="Only folders for specified K value will be created")
    args = parser.parse_args()
    return args


def generate_paths(base_path: str, folder_names: List[List[str]]) -> List[str]:
    """ Generates all combinations of folders specified in folder_names and prepends with *base_path*. Returns a list of generated paths. """
    name_combos = itertools.product(*folder_names)
    paths = list(map(lambda name_combo: path.join(base_path, *name_combo), name_combos))
    return paths


def generate_config_lines(folder_params: List[List[str]]) -> List[str]:
    """ Each folder has a specific config parameter associated with it. 
    In a similar to generate_paths fashion this function generates full combination of all params for each folder. """
    param_combos = itertools.product(*folder_params)
    config_lines = list(map(lambda param_combo: "\n".join(param_combo) + "\n", param_combos))
    return config_lines


def create_paths(target_folders: List[str]):
    """ Creates all paths specified in target_folders. """
    for folder in target_folders:
        target_path = Path(folder)
        if target_path.name == "basis" or target_path.name == "overlaps" or target_path.name == "eigencalc":
            target_path = target_path / ("out_" + target_path.name)
        if not path.exists(target_path):
            os.makedirs(target_path)


def multicopy_config(config_path: str, target_folders: List[str]) -> List[str]:
    """ Copies specified config file into specified list of directories. Returns a list of full paths to new configs. """
    config_name = path.basename(config_path)
    new_config_paths = list(map(lambda folder: path.join(folder, config_name), target_folders))
    for new_path in new_config_paths:
        shutil.copyfile(config_path, new_path)
    return new_config_paths


def set_config_params(config_paths: List[str], config_lines: List[str]):
    """ Appends i-th config_lines to i-th config in config_paths """
    for i in range(len(config_paths)):
        with open(config_paths[i], "a") as config:
            config.write(config_lines[i])


def main():
    args = parse_command_line_args()
    config_path = path.abspath(args.config)
    base_path = path.dirname(config_path)
    config = SpectrumSDTConfig(config_path)

    folder_names = [["K_" + args.K], ["even", "odd"], ["eigencalc", "properties"]]
    folder_params = [["K = " + args.K], ["symmetry = 0", "symmetry = 1"], ["stage = eigencalc", "stage = properties"]]

    if args.K.isdigit():
        # add extra stages if K is a single number
        folder_names[2].insert(0, "basis")
        folder_params[2].insert(0, "stage = basis")
        folder_names[2].insert(1, "overlaps")
        folder_params[2].insert(1, "stage = overlaps")
    else:
        # add parity if K is a range
        folder_names.insert(1, ["parity_0", "parity_1"])
        folder_params.insert(1, ["parity = 0", "parity = 1"])

    target_folders = generate_paths(base_path, folder_names)
    config_lines = generate_config_lines(folder_params)

    create_paths(target_folders)
    config_paths = multicopy_config(config_path, target_folders)
    set_config_params(config_paths, config_lines)


if __name__ == "__main__":
    main()
