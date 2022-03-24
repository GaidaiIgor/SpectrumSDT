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

    def get_mass_str(self) -> str:
        return self.params["mass"]

    def get_molecule_type(self) -> str:
        mass_string = self.get_mass_str()
        mass_tokens = mass_string.split(",")
        mass_tokens = [token.strip() for token in mass_tokens]
        # 2 types are supported
        if mass_tokens[0] == mass_tokens[2]:
            if mass_tokens[0] == mass_tokens[1]:
                return "AAA"
            else:
                return "ABA"
        else:
            raise Exception("Unknown molecule type")

    def get_J(self) -> int:
        return int(self.params["J"])

    def get_parity(self) -> int:
        return int(self.params["parity"]) if "parity" in self.params else None

    @staticmethod
    def parse_Ks(K_str, J, parity = None) -> List[int]:
        if K_str == "all":
            k_start = 1 if parity is None else (J + parity) % 2
            return [k_start, J]

        K_str_tokens = K_str.split("..")
        if len(K_str_tokens) == 1:
            return [int(K_str_tokens[0]), int(K_str_tokens[0])]
        else:
            return [int(K_str_tokens[0]), int(K_str_tokens[1])]

    def get_Ks(self) -> List[int]:
        return SpectrumSDTConfig.parse_Ks(self.params["K"], self.get_J(), self.get_parity())

    def get_symmetry(self) -> int:
        return int(self.params["basis"]["K0_symmetry"])

    def get_basis_fixed_J(self) -> int:
        return int(self.params["basis"]["fixed"]["J"])

    def get_basis_fixed_K(self) -> int:
        return int(self.params["basis"]["fixed"]["K"])

    def get_basis_fixed_root_path(self) -> str:
        return self.params["basis"]["fixed"]["root_path"]

    def get_number_of_states(self) -> int:
        mol_type = self.get_molecule_type()
        state_mult = 1 if mol_type == "AAA" else 3
        return int(self.params["eigensolve"]["num_states"]) * state_mult

    def get_grid_path(self) -> str:
        return self.params["grid_path"]

    def get_root_path(self) -> str:
        return self.params["root_path"]

    def get_geometric_phase(self) -> int:
        if "use_geometric_phase" not in self.params:
            return 0
        else:
            return int(self.params["use_geometric_phase"])

    def get_half_integers(self) -> bool:
        if "use_half_integers" not in self.params["basis"]:
            return False
        else:
            return bool(self.params["basis"]["use_half_integers"])

