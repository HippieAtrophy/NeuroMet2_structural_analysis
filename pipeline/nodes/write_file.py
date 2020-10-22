# -*- coding: utf-8 -*-


"""

"""
import os
import pandas as pd
import re
import shutil

from nipype import logging
from nipype.utils.filemanip import fname_presuffix, split_filename

from nipype import logging, LooseVersion
from nipype.utils.filemanip import fname_presuffix, check_depends
from nipype.interfaces.io import FreeSurferSource

from nipype.interfaces.base import (
    TraitedSpec,
    File,
    traits,
    Directory,
    InputMultiPath,
    OutputMultiPath,
    CommandLine,
    CommandLineInputSpec,
    isdefined,
    BaseInterfaceInputSpec,
    BaseInterface
)

from nipype.interfaces.freesurfer.base import FSCommand, FSTraitedSpec, FSTraitedSpecOpenMP, FSCommandOpenMP, Info
from nipype.interfaces.freesurfer.utils import copy2subjdir

__docformat__ = "restructuredtext"
iflogger = logging.getLogger("nipype.interface")

# Keeping this to avoid breaking external programs that depend on it, but
# this should not be used internally
FSVersion = Info.looseversion().vstring

class WriteFileInputSpec(BaseInterfaceInputSpec):

    data_str = traits.Str(argstr="%s", desc="subject id of surface file", mandatory=True)
    csv_file = File(desc="Output CSV file", mandatory=True)
    names = col_list = ['Measure:volume',
                         'Left-Lateral-Ventricle',
                         'Left-Inf-Lat-Vent',
                         'Left-Cerebellum-White-Matter',
                         'Left-Cerebellum-Cortex',
                         'Left-Thalamus',
                         'Left-Caudate',
                         'Left-Putamen',
                         'Left-Pallidum',
                         '3rd-Ventricle',
                         '4th-Ventricle',
                         'Brain-Stem',
                         'Left-Hippocampus',
                         'Left-Amygdala',
                         'CSF',
                         'Left-Accumbens-area',
                         'Left-VentralDC',
                         'Left-vessel',
                         'Left-choroid-plexus',
                         'Right-Lateral-Ventricle',
                         'Right-Inf-Lat-Vent',
                         'Right-Cerebellum-White-Matter',
                         'Right-Cerebellum-Cortex',
                         'Right-Thalamus',
                         'Right-Caudate',
                         'Right-Putamen',
                         'Right-Pallidum',
                         'Right-Hippocampus',
                         'Right-Amygdala',
                         'Right-Accumbens-area',
                         'Right-VentralDC',
                         'Right-vessel',
                         'Right-choroid-plexus',
                         '5th-Ventricle',
                         'WM-hypointensities',
                         'Left-WM-hypointensities',
                         'Right-WM-hypointensities',
                         'non-WM-hypointensities',
                         'Left-non-WM-hypointensities',
                         'Right-non-WM-hypointensities',
                         'Optic-Chiasm',
                         'CC_Posterior',
                         'CC_Mid_Posterior',
                         'CC_Central',
                         'CC_Mid_Anterior',
                         'CC_Anterior',
                         'BrainSegVol',
                         'BrainSegVolNotVent',
                         'VentricleChoroidVol',
                         'lhCortexVol',
                         'rhCortexVol',
                         'CortexVol',
                         'lhCerebralWhiteMatterVol',
                         'rhCerebralWhiteMatterVol',
                         'CerebralWhiteMatterVol',
                         'SubCortGrayVol',
                         'TotalGrayVol',
                         'SupraTentorialVol',
                         'SupraTentorialVolNotVent',
                         'MaskVol',
                         'BrainSegVol-to-eTIV',
                         'MaskVol-to-eTIV',
                         'lhSurfaceHoles',
                         'rhSurfaceHoles',
                         'SurfaceHoles',
                         'EstimatedTotalIntraCranialVol']

class WriteFileOutputSpec(TraitedSpec):

    csv_file = File(exists=True, desc="Output CSV file")


class WriteFile(BaseInterface):
    """
    ToDo: Example Usage:

    """

    input_spec =WriteFileInputSpec
    output_spec = WriteFileOutputSpec



    def _run_interface(self, runtime, correct_return_codes=(0,)):

        # if the csv file doesn't exists create one an write the title line:
        csv_file = self.inputs.csv_file
        if not os.path.isfile(csv_file):
            with open(csv_file, 'w+') as f:
                f.write(', '.join(self.inputs.names) + ' \n')
        # then write the data string.
        with open(csv_file, 'a') as f:
            f.write(self.inputs.data_str + ' \n')
        return runtime

    def _list_outputs(self):

        outputs = self._outputs().get()
        outputs["csv_file"] = self.inputs.csv_file
        return outputs
