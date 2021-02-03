import os.path as path
from typing import Any, Dict, List, Tuple

class SpectrumSDTConfig:
    """ Provides access to information in config file """

    def read_inner_dict(self, config_lines, line_num) -> Tuple[Dict[str, Any], int]:
        dict = {}
        while line_num + 1 < len(config_lines):
            line_num += 1
            line = config_lines[line_num]
            line = line.split("!")[0].strip() # remove comments
            if len(line) == 0:
                continue
            if line == ")":
                return dict, line_num
            tokens = line.split("=")
            key = tokens[0].strip()
            value = tokens[1].strip()
            if value == "(":
                value, line_num = self.read_inner_dict(config_lines, line_num)
            dict[key] = value
        return dict, line_num

    def __init__(self, config_path):
        self.config_path = config_path
        with open(config_path) as config_file:
            config_lines = config_file.readlines()
        self.params, _ = self.read_inner_dict(config_lines, -1)

    def get_stage(self) -> str:
        return self.params["stage"]

    def get_fixed_basis_JK(self) -> int:
        if "use_fixed_basis_JK" in self.params:
            return int(self.params["use_fixed_basis_JK"])
        else:
            return 0

    def get_mass_str(self) -> str:
        return self.params["mass"]

    def get_J(self) -> int:
        return int(self.params["J"])

    def get_Ks(self) -> List[int]:
        K_str = self.params["K"]
        if K_str == "all":
            J = self.get_J()
            k_start = (J + self.get_parity()) % 2
            return [k_start, J]

        K_str_tokens = K_str.split("..")
        if len(K_str_tokens) == 1:
            return [int(K_str_tokens[0]), int(K_str_tokens[0])]
        else:
            return [int(K_str_tokens[0]), int(K_str_tokens[1])]

    def get_parity(self) -> int:
        return int(self.params["parity"])

    def get_symmetry(self) -> int:
        return int(self.params["basis"]["symmetry"])

    def get_basis_root_path(self) -> str:
        return self.params["basis_root_path"]

    def get_basis_J(self) -> int:
        return int(self.params["basis_J"])

    def get_basis_K(self) -> int:
        return int(self.params["basis_K"])

    def get_number_of_states(self) -> int:
        return int(self.params["eigencalc"]["num_states"])

    def get_ncv(self) -> int:
        return int(self.params["ncv"])

    def get_grid_path(self) -> str:
        return self.params["grid_path"]

    def get_root_path(self) -> str:
        return self.params["root_path"]
