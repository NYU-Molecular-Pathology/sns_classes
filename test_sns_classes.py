#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
unit tests for the find module
"""
import sys
import unittest
import os
import yaml
import csv
from collections import defaultdict
from classes import AnalysisItemMissing
from classes import AnalysisItem
from classes import SnsWESAnalysisOutput
from classes import SnsAnalysisSample
import config

# add parent dir to sys.path to import util
sys.path.insert(0, "..")
from util import log
sys.path.pop(0)

scriptdir = os.path.dirname(os.path.realpath(__file__))
fixture_dir = os.path.join(scriptdir, "fixtures")
sns_output_dir = os.path.join(fixture_dir, 'sns_output')
sns_analysis1_dir = os.path.join(sns_output_dir, 'sns_analysis1')
sns_analysis1_nosettings_dir = os.path.join(sns_output_dir, 'sns_analysis1_nosettings')
sns_analysis1_qsuberrors_dir = os.path.join(sns_output_dir, "sns_analysis1_qsuberrors")
sns_analysis1_summaryX_dir = os.path.join(sns_output_dir, "sns_analysis1_summaryX")


# ~~~~ SETUP CONFIGS FROM EXTERNAL RESOURCES ~~~~~~ #
config_dir = "config"
config_file = os.path.join(scriptdir, config_dir, "sns.yml")
with open(config_file, "r") as f:
    configs = yaml.load(f)



class TestAnalysisItem(unittest.TestCase):
    def setUp(self):
        self.analysis_item1 = AnalysisItem(id = 'foo')
        # remove the handlers because they are too verbose here
        self.analysis_item1.logger = log.remove_all_handlers(logger = self.analysis_item1.logger)

    def tearDown(self):
        del self.analysis_item1

    def test_true(self):
        self.assertTrue(True, 'Demo True assertion')

    def test_AnalysisItem_files_type(self):
        """
        Test a NextSeq demo run that should be valid
        """
        x = self.analysis_item1
        self.assertTrue(isinstance(x.files, defaultdict), "AnalysisItem 'files' is not type 'defaultdict'")

    def test_AnalysisItem_files_entry1(self):
        """
        Test a NextSeq demo run that should be valid
        """
        x = self.analysis_item1
        self.assertTrue(isinstance(x.files[0], list), "AnalysisItem 'files' entry is not type 'list'")

    def test_AnalysisItem_files_entry_missingkey(self):
        """
        Test a NextSeq demo run that should be valid
        """
        item = self.analysis_item1.get_files('doesntexist')
        self.assertTrue(isinstance(item, list), "AnalysisItem 'get_files' method key did not return type 'list' for missing key")

    def test_AnalysisItem_files_entry_listnone(self):
        """
        Test a NextSeq demo run that should be valid
        """
        item = self.analysis_item1.list_none(l = self.analysis_item1.get_files('doesntexist'))
        self.assertIsNone(item, "AnalysisItem 'list_none' method key did not return type 'None' for empty list")


class TestSnsWESAnalysisOutput(unittest.TestCase):
    def setUp(self):
        self.analysis_output_1 = SnsWESAnalysisOutput(dir = sns_analysis1_dir, id = 'sns_analysis1', sns_config = configs)
        # remove the handlers because they are too verbose here
        self.analysis_output_1.logger = log.remove_all_handlers(logger = self.analysis_output_1.logger)

        self.sns_analysis1_nosettings = SnsWESAnalysisOutput(dir = sns_analysis1_nosettings_dir, id = 'analysis1_nosettings', sns_config = configs, debug = True)
        self.sns_analysis1_nosettings.logger = log.remove_all_handlers(logger = self.sns_analysis1_nosettings.logger)

        self.sns_analysis1_qsuberrors = SnsWESAnalysisOutput(dir = sns_analysis1_qsuberrors_dir, id = 'sns_analysis1_qsuberrors', sns_config = configs, debug = True)
        self.sns_analysis1_qsuberrors.logger = log.remove_all_handlers(logger = self.sns_analysis1_qsuberrors.logger)

        self.sns_analysis1_summaryX = SnsWESAnalysisOutput(dir = sns_analysis1_summaryX_dir, id = 'sns_analysis1_summaryX', sns_config = configs, debug = True)
        self.sns_analysis1_summaryX.logger = log.remove_all_handlers(logger = self.sns_analysis1_summaryX.logger)


    def tearDown(self):
        del self.analysis_output_1
        del self.sns_analysis1_nosettings
        del self.sns_analysis1_qsuberrors
        del self.sns_analysis1_summaryX

    def test_invalid_path(self):
        analysis_dir = os.path.join(sns_output_dir, 'foo')
        # x = SnsWESAnalysisOutput(dir = analysis_dir, id = 'foo', sns_config = configs)
        # x.logger = log.remove_all_handlers(logger = x.logger)
        # self.assertFalse(x.validate(), 'Analysis dir with invalid path returns True validation')
        with self.assertRaises(AnalysisItemMissing):
            SnsWESAnalysisOutput(dir = analysis_dir, id = 'foo', sns_config = configs)

    def test_no_settings(self):
        self.assertRaises(AnalysisItemMissing, self.sns_analysis1_nosettings.validate)


    def test_qsub_errors(self):
        self.assertFalse(self.sns_analysis1_qsuberrors.validate(), 'Analysis dir with no settings file returned True validation')

    def test_summary_combined_errors(self):
        summary_combined_wes_file = self.sns_analysis1_summaryX.static_files.get('summary_combined_wes', None)
        rows = []
        with open(summary_combined_wes_file, 'rb') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                rows.append(row)
        contains_errors = self.sns_analysis1_summaryX.summary_combined_contains_errors(summary_combined_wes_rows = rows, err_pattern = 'X')
        self.assertTrue(contains_errors)

    def test_valid_analysis_output(self):
        self.assertTrue(self.analysis_output_1.validate(), 'Valid analysis dir returned False validation')

    def test_get_samples(self):
        self.assertTrue(len(self.analysis_output_1.get_samples()) == 4, 'Analysis dir "analysis_output_1" did not return 4 samples')

    def test_get_bam_dir_exists(self):
        self.assertTrue(bool(self.analysis_output_1.dirs.get('BAM-DD', None)), 'Analysis dir "analysis_output_1" did not return an entry for "BAM-DD" samples')



if __name__ == '__main__':
    unittest.main()
