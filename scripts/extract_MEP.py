#!/usr/bin/env python3

""" This script extracts minimum energy path (MEP) along rho from given potential and grids files.
Launch this script from within the folder with grid files and pes_out.txt.
The result will be written to a file named MEP_rho.txt.
This file can be used to generate optimized rho-grid (see manual for details). """


def main():
    grid_sizes = [0]*3
    with open("grid_rho.txt") as grid_file:
        grid_sizes[0] = int(grid_file.readline().split()[3])
        grid_rho = [float(line.split()[0]) for line in grid_file]
    with open("grid_theta.txt") as grid_file:
        grid_sizes[1] = int(grid_file.readline().split()[3])
    with open("grid_phi.txt") as grid_file:
        grid_sizes[2] = int(grid_file.readline().split()[3])
    with open("pes_out.txt") as potential_file:
        potential = [float(line) for line in potential_file]
    rho_step = grid_sizes[1]*grid_sizes[2]
    mep_rho = [min(potential[rho_ind*rho_step : (rho_ind+1)*rho_step]) for rho_ind in range(grid_sizes[0])]
    with open("MEP_rho.txt", "w") as mep_file:
        for rho_ind in range(grid_sizes[0]):
            mep_file.write("{0} {1}\n".format(grid_rho[rho_ind], mep_rho[rho_ind]))


if __name__ == "__main__":
    main()
