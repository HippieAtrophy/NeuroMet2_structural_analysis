# -*- coding: utf-8 -*-


"""
Emulate Freesurfer QDec wrapping *stats2table commands
"""
import os
import glob
import subprocess
from nipype import logging
from nipype.interfaces.base import (
    traits,
    TraitedSpec,
    Directory,
    BaseInterfaceInputSpec,
    BaseInterface
)
from traits.trait_types import Str

__docformat__ = "restructuredtext"
iflogger = logging.getLogger("nipype.interface")


class QDecInputSpec(BaseInterfaceInputSpec):

    basedir = traits.Str(desc="Base directory", argstr="%s", mandatory=True)
    fs_dir_template = traits.Str("*/*.freesurfer", argstr="%s", desc="Freesurfer directory template, default sub*/*.freesurfer", usedefault = True)


class QDecOutputSpec(TraitedSpec):

    stats_directory = Directory(desc="stat_tables directory", mandatory=True)
    #outputs = traits.List(traits.Str(), desc="stdout messages")
    #errors = traits.List(traits.Str(), desc="stderr messages")


class QDec(BaseInterface):
    """
    ToDo: Example Usage:

    """

    input_spec = QDecInputSpec
    output_spec = QDecOutputSpec

    def __make_sublist(self):
        sublist = glob.glob(self.input_spec.basedir + '/*/')
        return sublist

    def _run_interface(self, runtime, correct_return_codes=(0,)):

        os.environ["SUBJECTS_DIR"] = str(self.input_spec.basedir)
        outputs = list()
        errors = list()
        commands =["asegstats2table --common-segs --meas volume --tablefile {d}/aseg.volume.stats.dat --statsfile=aseg.stats --subjects {s}",
                   "asegstats2table --common-segs --meas volume --tablefile {d}/wmparc.volume.stats.dat --statsfile=wmparc.stats --subjects {s}",
                   "aparcstats2table --hemi lh --parc aparc --meas area --tablefile {d}/lh.aparc.area.stats.dat --subjects {s}",
                   "aparcstats2table --hemi lh --parc aparc --meas volume --tablefile {d}/lh.aparc.volume.stats.dat --subjects {s}",
                   "aparcstats2table --hemi lh --parc aparc --meas thickness --tablefile {d}/lh.aparc.thickness.stats.dat --subjects {s}",
                   "aparcstats2table --hemi lh --parc aparc --meas meancurv --tablefile {d}/lh.aparc.meancurv.stats.dat --subjects {s}",
                   "aparcstats2table --hemi lh --parc aparc.a2009s --meas area --tablefile {d}/lh.aparc.a2009s.area.stats.dat --subjects {s}",
                   "aparcstats2table --hemi lh --parc aparc.a2009s --meas volume --tablefile {d}/lh.aparc.a2009s.volume.stats.dat --subjects {s}",
                   "aparcstats2table --hemi lh --parc aparc.a2009s --meas thickness --tablefile {d}/lh.aparc.a2009s.thickness.stats.dat --subjects {s}",
                   "aparcstats2table --hemi lh --parc aparc.a2009s --meas meancurv --tablefile {d}/lh.aparc.a2009s.meancurv.stats.dat --subjects {s}",
                   "aparcstats2table --hemi rh --parc aparc --meas area --tablefile {d}/rh.aparc.area.stats.dat --subjects {s}",
                   "aparcstats2table --hemi rh --parc aparc --meas volume --tablefile {d}/rh.aparc.volume.stats.dat --subjects {s}",
                   "aparcstats2table --hemi rh --parc aparc --meas thickness --tablefile {d}/rh.aparc.thickness.stats.dat --subjects {s}",
                   "aparcstats2table --hemi rh --parc aparc --meas meancurv --tablefile {d}/rh.aparc.meancurv.stats.dat --subjects {s}",
                   "aparcstats2table --hemi rh --parc aparc.a2009s --meas area --tablefile {d}/rh.aparc.a2009s.area.stats.dat --subjects {s}",
                   "aparcstats2table --hemi rh --parc aparc.a2009s --meas volume --tablefile {d}/rh.aparc.a2009s.volume.stats.dat --subjects {s}",
                   "aparcstats2table --hemi rh --parc aparc.a2009s --meas thickness --tablefile {d}/rh.aparc.a2009s.thickness.stats.dat --subjects {s}",
                   "aparcstats2table --hemi rh --parc aparc.a2009s --meas meancurv --tablefile {d}/rh.aparc.a2009s.meancurv.stats.dat --subjects {s}"
                   ]
        for command in commands:
            process = subprocess.Popen(command.format(d=os.path.join(str(self.input_spec.basedir), 'stats_tables'),
                                                      s=' '.join(self.__make_sublist())).split(),
                                        stdout=subprocess.PIPE)
            output, error = process.communicate()
            outputs.append(output)
            errors.append(error)
        setattr(self, '_outputs', outputs)  # Save result
        setattr(self, '_errors', errors)  # Save result

        return runtime

    def _list_outputs(self):

        outputs = self._outputs().get()
        outputs["stats_directory"] = os.path.join(str(self.input_spec.basedir), 'stats_tables')
        #outputs["outputs"] = getattr(self, '_outputs')
        #outputs["errors"] = getattr(self, '_errors')
        return outputs
