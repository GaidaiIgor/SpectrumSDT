#!/usr/bin/env python3

import argparse
import itertools
import os
import shutil
import os.path as path
from pathlib import Path
from typing import List, Dict
from SpectrumSDTConfig import SpectrumSDTConfig


def parse_command_line_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Creates a folder structure for a given value of K, copies config template to the target folders, and fills out placeholder values. "
            "Placeholder values can be referred to in config template as {parameter}, where parameter name can be one of: K, symmetry, stage or parity (if K is a range).")
    parser.add_argument("-c", "--config", default="spectrumsdt.config", help="Path to a template configuration file. ./spectrumsdt.config by default.")
    parser.add_argument("-K", required=True, help="The value of K for which to generate the folder structure. "
        "Can either be a scalar, a custom range in the form K1..K2, or 'all' to include all Ks for a given J and parity.")
    parser.add_argument("--extra", help="A dictionary with extra placeholder values that are applied to all target folders.")
    args = parser.parse_args()
    return args


def resolve_defaults(args: argparse.Namespace):
    if args.extra is None:
        args.extra = {}
    else:
        args.extra = eval(args.extra)
        for key, value in args.extra.items():
            if not isinstance(value, str):
                args.extra[key] = str(value)


def get_available_symmetries(molecule_type: str, use_geometric_phase: int) -> List[str]:
    if (molecule_type == "AAA"):
        if (use_geometric_phase == 0):
            return ["0", "1", "2", "3"]
        else:
            return ["0", "2"]
    else:
        if (use_geometric_phase == 0):
            return ["0", "1"]
        else:
            return ["0"]


def generate_paths(base_path: str, param_names: List[str], use_param_names: List[bool], value_combos) -> List[str]:
    """ Generates all config paths for all combination of values of all parameters. """
    all_paths = []
    for i in range(len(value_combos)):
        next_path = base_path
        next_combo = value_combos[i]
        for j in range(len(next_combo)):
            if use_param_names[j]:
                next_path = path.join(next_path, param_names[j] + "_" + next_combo[j])
            else:
                next_path = path.join(next_path, next_combo[j])
        all_paths.append(next_path)
    return all_paths


def create_paths(target_paths: List[str]):
    """ Creates all paths specified in target_paths. """
    for next_path in target_paths:
        next_path = Path(next_path)
        create_bin = any(list(map(lambda part: part == "basis" or part == "overlaps" or part == "eigensolve", next_path.parts[-2:])))
        if create_bin:
            next_path = next_path / "bin"
        if not path.exists(next_path):
            os.makedirs(next_path)


def multicopy_config(config_path: str, target_paths: List[str]) -> List[str]:
    """ Copies specified config file into specified list of directories. """
    for target_path in target_paths:
        shutil.copyfile(config_path, path.join(target_path, "spectrumsdt.config"))


def generate_param_dicts(param_names: List[str], value_combos, extra: Dict[str, str]) -> List[Dict[str, str]]:
    """ Transforms lists of names and value combinations to list of dictionaries. """
    res = []
    for i in range(len(value_combos)):
        next_dict = {}
        for j in range(len(value_combos[i])):
            next_dict[param_names[j]] = value_combos[i][j]
        next_dict.update(extra)
        if not 'parity' in next_dict:
            # If parity is not defined, then it does not matter, set it to 0
            next_dict['parity'] = 0
        res.append(next_dict)
    return res


def set_placeholder_params(target_paths: List[str], param_dicts: List[Dict[str, str]]):
    """ Iterates over configs at *target_paths* and sets placeholder params according to *param_dicts*. """
    for i in range(len(target_paths)):
        with open(path.join(target_paths[i], "spectrumsdt.config"), "r+") as config:
            content = config.read()
            formatted = content.format(**param_dicts[i])
            config.seek(0)
            config.write(formatted)
            config.truncate()


def main():
    args = parse_command_line_args()
    resolve_defaults(args)
    config_path = path.abspath(args.config)
    config = SpectrumSDTConfig(config_path)

    molecule_type = config.get_molecule_type()
    use_geometric_phase = config.get_geometric_phase()
    available_symmetries = get_available_symmetries(molecule_type, use_geometric_phase)

    # this is general template, which may be modified below
    param_names = ["K", "symmetry", "stage"]
    use_param_names = [True, True, False]
    param_values = [[args.K], available_symmetries, ["eigensolve", "properties"]]

    Ks = SpectrumSDTConfig.parse_Ks(args.K, config.get_J())
    if Ks[0] == Ks[1]:
        do_basis = False
        if "fixed" in config.params["basis"]:
            J = config.get_J()
            fixed_J = config.get_basis_fixed_J()
            fixed_K = config.get_basis_fixed_K()
            if J == fixed_J and Ks[0] == fixed_K:
                do_basis = True
        else:
            do_basis = True

        if do_basis:
            # add extra stages if basis calculation is needed
            param_values[2].insert(0, "basis")
            param_values[2].insert(1, "overlaps")

    if config.params["use_rovib_coupling"] == "1" and Ks[0] <= 1:
        param_names.insert(3, "parity")
        use_param_names.insert(3, True)
        if Ks[0] == 0:
            param_values.insert(3, str(config.get_J() % 2))
        else:
            param_values.insert(3, ["0", "1"])

    if len(param_values) == 4 and len(param_values[2]) == 4:
        parity_insensitive_values = param_values[0:2]
        parity_insensitive_values[2] = parity_insensitive_values[2][0:1]
        parity_sensitive_values = param_values
        parity_sensitive_values[2] = parity_sensitive_values[2][2:3]
        value_combos = list(itertools.product(*parity_insensitive_values)) + list(itertools.product(*parity_sensitive_values))
    else:
        value_combos = list(itertools.product(*param_values))

    base_path = path.dirname(config_path)
    target_paths = generate_paths(base_path, param_names, use_param_names, value_combos)
    create_paths(target_paths)
    multicopy_config(config_path, target_paths)

    param_dicts = generate_param_dicts(param_names, value_combos, args.extra)
    set_placeholder_params(target_paths, param_dicts)


if __name__ == "__main__":
    main()
