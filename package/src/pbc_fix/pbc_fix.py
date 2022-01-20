VERSION = '2.0.6'

import argparse
import MDAnalysis
import sys

import numpy as np
import numpy.typing as npt
from pathlib import Path, PosixPath
from typing import Optional

from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.core.universe import Universe
from MDAnalysis.coordinates.base import Timestep
from MDAnalysis.coordinates.XTC import XTCWriter
from MDAnalysis.transformations import center_in_box

#NOTE: if executing pbc_fix.py directly the import statement below
# must be modified to remove the relative import:
# from utils.trajectory_processsing import TrajectorySparseIterator
from .utils.trajectory_processing import TrajectorySparseIterator


class PBCFixOutputs:
    def __init__(
        self,
        output_dir: PosixPath,
        output_name: str,
        num_atoms: int,
        t_start: int,
        gro_output_frequency: int,
    ):
        # Create output dir
        output_dir.mkdir(parents=True, exist_ok=False)
        print(f'Writing outputs at {output_dir}')

        # Add log file containing command
        with open(output_dir.joinpath('pbc_fix.log'), 'w') as f:
            f.write(f'[pbc_fix {VERSION}]\n')
            f.write('\nThis folder and its content were created using the following command:\n\n')
            f.write(' '.join(sys.argv) + '\n')

        # Create output path stem for output files
        self.output_name = output_dir.joinpath(f'{output_name}_fixed')

        # Create xtc writer
        self.traj_writer = XTCWriter(str(self.output_name.with_suffix('.xtc')), num_atoms)

        # Store parameters for outputting
        self.t_start = t_start
        self.gro_output_frequency = gro_output_frequency
        self.gro_outputs = 0

    def get_gro_output(self, time: int, force: bool = False):
        if force or (time - self.t_start) // self.gro_output_frequency > self.gro_outputs:
            self.gro_outputs += 1
            return str(self.output_name.with_name(f'{self.output_name.stem}_ref{int(time/1000)}ns.gro'))
        return None


def parse_args():
    parser = argparse.ArgumentParser(
        prog='pbc_fix',
        usage='see pbc_fix --help',
        description='Typical usage: pbc_fix --topology top.gro --trajectory traj.xtc',
        add_help=True,
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=28, width=100)
    )

    # Input/output
    parser.add_argument(
        '--trajectory', required=True, metavar='[required]', help='Trajectory file (.xtc) for which to fix pbc.'
    )
    parser.add_argument(
        '--topology',
        required=True,
        metavar='[required]',
        help='Topology file (e.g. .gro) used to load the trajectory.',
    )
    parser.add_argument(
        '--output',
        metavar='[optional]',
        help='Directory where to save the output of the script. Defaults to the name of the trajectory file.',
    )
    parser.add_argument(
        '--gro-output-frequency',
        default=-1,
        type=int,
        metavar='[optional]',
        help='If specified, time interval (in ns) at which to output a gro file. Useful to check the script is' \
            + ' working as expected. The gro file corresponding to the last processed frame is *always* generated.',
    )

    # Selection of particles of interest
    parser.add_argument(
        '--selection',
        default='not resname W and not resname ION',
        metavar='[optional]',
        help='Selection string identifying the particles of interest in the trajectory. The outputs produced will' \
            + ' *only* contain this selection. Defaults to \"not resname W and not resname ION\".'
    )

    # Reference structure
    parser.add_argument(
        '--ref-topology',
        action='store_true',
        help='By default the reference coordinates of the selection are taken from the first frame of the trajectory.' \
            + ' Use this flag to take them from the topology file instead.',
    )

    # Coordinates modification parameters
    parser.add_argument(
        '--buffer',
        default=80,
        type=float,
        metavar='[optional]',
        help='Buffer (in Å) to maintain between the lowest/highest particules and the periodic box bottom/top.' \
            + 'Defaults to 80.',
    )

    # Version info
    parser.add_argument('--version', action='version', version=f'pbc_fix {VERSION}')

    return parser.parse_args()


def process_args(args):
    args.gro_output_frequency *= 1000

    # Check buffer is positive
    if args.buffer < 0:
        raise ValueError(f'--buffer should be positive. Got {args.buffer}')

    # Check output dir
    args.topology = Path(args.topology)
    args.trajectory = Path(args.trajectory)
    if args.output is None:
        args.output = args.trajectory.parent.joinpath(args.trajectory.stem)
    args.output = Path(args.output)
    if args.output.exists():
        raise ValueError(f'Output directory {args.output} already exists.')

    return args


def get_reference_coordinates(
    sel: str,
    u: Optional[Universe] = None,
    topology: Optional[str] = None,
    ref_topology: bool = False,
):
    if ref_topology:
        ag = MDAnalysis.Universe(topology).select_atoms(sel)
    else:
        u.trajectory[0]
        ag = u.select_atoms(sel)
    return ag.positions[:, 2]


def update_frame(
    ts: Timestep,
    ag: AtomGroup,
    z_ref: npt.ArrayLike,
    z_buffer: float,
    output_traj: Optional[str] = None,
    output_gro: Optional[str] = None,
):
    # Get current box dimension along z axis
    box_dim = ts.dimensions[2]

    # Get current coordinates
    coord = ag.positions

    # Calculate displacement
    coord_z_delta = coord[:, 2] - z_ref

    # This works because we'll only update the coord once the delta is greater than half the box
    # NOTE: this is a vector
    n_box_z = (np.abs(coord_z_delta) + 0.5 * box_dim) // box_dim

    # Update coords: deal with pbc
    to_move_up = coord_z_delta < 0  # these particles are candidates for crossing the upper boundary of the box
    to_move_down = coord_z_delta > 0  # these particles are candidates for crossing the lower boundary of the box
    coord[to_move_up, 2] += n_box_z[to_move_up] * box_dim
    coord[to_move_down, 2] -= n_box_z[to_move_down] * box_dim

    # Store coordinates for next frame comparison
    z_ref = coord[:, 2]

    # Update box size: encompass all coords and maintain buffer
    box_buffer = box_dim - (np.max(coord[:, 2]) - np.min(coord[:, 2]))
    if box_buffer < z_buffer:
        box_dim += z_buffer - box_buffer

    # Apply changes
    ts.dimensions[2] = box_dim
    ag.positions = coord

    # Center coords in box
    tf = center_in_box(ag)
    tf(ts)

    # write frame
    if output_gro is not None:
        ag.write(output_gro)

    # write trajectory
    if output_traj is not None:
        output_traj.write(ag)

    return z_ref


def process_trajectory(
    path_topology: PosixPath,
    path_trajectory: PosixPath,
    path_out: PosixPath,
    sel: str,
    pbc_buffer: int,
    gro_output_frequency: int,
    ref_topology: bool = False,
):
    # Load trajectory
    u = MDAnalysis.Universe(str(path_topology), str(path_trajectory))

    # Select particles of interest
    ag = u.select_atoms(sel)

    # Get initial reference coordinates
    z_ref = get_reference_coordinates(sel, u, path_topology, ref_topology)

    # Create output writer
    pbc_outputs = PBCFixOutputs(
        output_dir=path_out,
        output_name=path_trajectory.stem,
        num_atoms=ag.n_atoms,
        t_start=u.trajectory[0].time,
        gro_output_frequency=gro_output_frequency,
    )

    # Process frames
    traj_iterator = TrajectorySparseIterator(u.trajectory)
    for ts in traj_iterator:
        z_ref = update_frame(
            ts=ts,
            ag=ag,
            z_ref=z_ref,
            z_buffer=pbc_buffer,
            output_traj=pbc_outputs.traj_writer,
            output_gro=pbc_outputs.get_gro_output(time=ts.time, force=traj_iterator.last_frame_reached),
        )

    print(f'Done. Check results in {pbc_outputs.output_name}.')


def main():
    args = parse_args()
    args = process_args(args)
    process_trajectory(
        path_topology=args.topology,
        path_trajectory=args.trajectory,
        path_out=args.output,
        sel=args.selection,
        pbc_buffer=args.buffer,
        gro_output_frequency=args.gro_output_frequency,
        ref_topology=args.ref_topology,
    )


if __name__ == '__main__':
    main()
