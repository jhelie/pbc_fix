import sys
import numpy as np
from typing import Optional
from MDAnalysis.coordinates.XTC import XTCReader

class TrajectorySparseIterator:
    def __init__(
        self,
        traj: XTCReader,
        frame_step: int = 1,
        t_start: Optional[int] = None,
        t_end: Optional[int] = None,
        verbose=True,
    ):
        self.traj = traj
        self.f_start = 0 if t_start is None else int((t_start - traj[0].time) / traj.dt)
        self.f_end = traj.n_frames - 1 if t_end is None else int((t_end - traj[0].time) / traj.dt)
        self.f_step = frame_step
        self.frames_to_process = np.arange(self.f_start, self.f_end, self.f_step)
        self.num_frames_to_process = len(self.frames_to_process)
        self.num_frames_trajectory = self.traj.n_frames
        self.verbose = verbose
        self.counter = 0
        self.last_frame_reached = self.counter == self.num_frames_to_process

    def __iter__(self):
        # Re-initialise variables used for iteration so that the iterator can be re-used
        self.counter = 0
        self.last_frame_reached = self.counter == self.num_frames_to_process
        return self

    def __next__(self):
        if self.last_frame_reached:
            raise StopIteration

        # Update counter
        self.counter += 1

        # Update boolean indicating whether the last frame has been reached
        self.last_frame_reached = self.counter == self.num_frames_to_process

        # Display message about progress of trajectory parsing
        if self.verbose:
            msg = f'\rProcessing frame {self.counter}/{self.num_frames_to_process}'
            msg += f' (every {self.f_step} from {self.f_start} to {self.f_end} out of {self.num_frames_trajectory})'
            sys.stdout.flush()
            sys.stdout.write(msg)

        # Return TimeStep corresponding to frame of interest
        return self.traj[self.frames_to_process[self.counter - 1]]
