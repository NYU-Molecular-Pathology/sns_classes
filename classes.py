#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Object classes for abstracting & interacting with ``sns`` pipeline output
"""
import os
import sys
import csv
import json
from collections import defaultdict

# add parent dir to sys.path to import util
scriptdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(scriptdir)
sys.path.insert(0, parentdir)
from util import find
from util import log
from util import tools
from util.classes import LoggedObject
from util.classes import AnalysisItem
sys.path.pop(0)
import config # config local to this dir

# ~~~~ CUSTOM CLASSES ~~~~~~ #
class AnalysisItemMissing(Exception):
    """
    An exception to use if an item from the analysis is missing

    Examples
    --------
    Example usage::

        raise AnalysisItemMissing(message = err_message, errors = '')

    """
    def __init__(self, message, errors):
        super(AnalysisItemMissing, self).__init__(message)
        self.errors = errors

class AnalysisInvalid(Exception):
    """
    An exception to use if the analysis is determined to be invalid

    Examples
    --------
    Example usage::

        raise AnalysisInvalid(message = err_message, errors = '')

    """
    def __init__(self, message, errors):
        super(AnalysisInvalid, self).__init__(message)
        self.errors = errors


class SnsWESAnalysisOutput(AnalysisItem):
    """
    Container for metadata about a sns WES targeted exome sequencing run analysis
    """
    def __init__(self, dir, id, sns_config = None, results_id = None, extra_handlers = None, debug = False):
        """
        Parameters
        ----------
        dir: str
            path to the directory containing ``sns`` pipeline output
        id: str
            ID for the analysis, typically the parent analysis output dir name, corresponding to a NextSeq run ID
        results_id: str
            typically a time-stamped ID of the results for the analysis, and the subdir name for the anaysis output
        sns_config: dict
            configuration items for the run; requires 'analysis_output_index' dict, and 'email_recipients'
        extra_handlers: list
            a list of Filehandlers, or ``None``
        debug: bool
            whether the analysis output should be intitialized in `debug` mode which skips validation

        Examples
        --------
        Example usage::

            import config
            from sns_classes import SnsWESAnalysisOutput
            analysis_dir = "/ifs/data/molecpathlab/NGS580_WES/170623_NB501073_0015_AHY5Y3BGX2/results_2017-06-26_20-11-26"
            analysis_id = "170623_NB501073_0015_AHY5Y3BGX2"
            results_id = "results_2017-06-26_20-11-26"
            analysis = SnsWESAnalysisOutput(dir = analysis_dir, id = analysis_id, results_id = results_id, sns_config = config.sns)

            import config
            from sns_classes import SnsWESAnalysisOutput
            analysis_dir = '/ifs/data/molecpathlab/scripts/snsxt/snsxt/fixtures/sns_output/sns_analysis1'
            analysis_id = "sample_analysis"
            results_id = "results1"
            analysis = SnsWESAnalysisOutput(dir = analysis_dir, id = analysis_id, results_id = results_id, sns_config = config.sns)

        Todo
        ----
        Reverted ``sns_config`` to being loaded from the config local to this file's dir. Need to deprecate this config entirely. Do not use it anymore.


        """
        AnalysisItem.__init__(self, id = id, extra_handlers = extra_handlers)
        # ID for the analysis run output; should match NextSeq ID
        self.id = str(id)
        self._init_logger(extra_handlers = extra_handlers)

        # path to the directory containing analysis output
        self.dir = os.path.abspath(dir)

        # config dict for sns program settings
        if not sns_config:
            # if not config passed, used default module
            self.sns_config = config.sns
        else:
            self.sns_config = sns_config
        # check that 'analysis_output_index' is in sns_config, else use default one
        if not self.sns_config.get('analysis_output_index', None):
            self.logger.debug('"analysis_output_index" key not found in passed "sns_config" dict, loading default sns_config instead')
            self.sns_config = config.sns
        # TODO: this is deprecated, need to remove it! Dont use it!!


        # timestamped ID for the analysis results, if supplied
        self.results_id = str(results_id)

        self._init_attrs()
        self._init_dirs()
        self._init_files()
        self._init_static_files()

        # self._init_analysis_config()

        # get the samples for the analysis
        # self.samples = self.get_samples()

        # the object should try to validate itself upon initialization
        # validation will fail if some static files are not present
        # this should kill the program, since it means the analysis output is invalid
        # maybe change this later if needed
        # try:
        if not debug: # True = dont validate, False = validate
            self.is_valid = self.validate()
            # if not self.is_valid:
            #     err_message = 'Analysis did not pass validations:\n{0}'.format(self.__repr__())
            #     raise AnalysisInvalid(message = err_message, errors = '')
            # # Don't raise here, because it already gets raised in run.py



    def __repr__(self):
        return("SnsWESAnalysisOutput {0} ({1}) located at {2}".format(self.id, self.results_id, self.dir))

    def _init_logger(self, extra_handlers = None):
        """
        Initializes the logger for the object

        Todo
        ----
        This appears to be causing double log messages, review this logger and determine if its necessary to keep it
        """
        # extra log handlers
        self.extra_handlers = extra_handlers
        # set up per-analysis logger
        self.logger = log.build_logger(name = self.id)
        if self.extra_handlers:
            self.logger = log.add_handlers(logger = self.logger, handlers = extra_handlers)
        self.logger.debug("Initialized logging for analysis: {0}".format(self.id))

    def _init_attrs(self):
        """
        Initializes attributes for the analysis

        Todo
        ----
        This is obtaining configs from the local config file; don't use these configs anymore, dont use these attributes, need to remove them
        """
        self.email_recipients = self.sns_config['email_recipients']
        self.analysis_output_index = self.sns_config['analysis_output_index']

    def _init_dirs(self):
        """
        Initializes the path attributes for items associated with the sequencing run
        from list of dirnames and filename patterns for the output steps in the sns WES analysis output

        Todo
        ----
        This is obtaining configs from the local config file; don't use these configs anymore, dont use these attributes, need to remove them. When using this module with ``snsxt``, the tasks should instead get the files explicitly from the tasks' ``input_dir``
        """
        for name, attributes in self.analysis_output_index.items():
            if name not in ['_parent']:
                self.set_dir(name = name, path = find.find(search_dir = self.dir, inclusion_patterns = name, search_type = "dir", num_limit = 1, level_limit = 0))

    def _init_static_files(self):
        """
        Initializes paths to files that should always exist in the same location for an analysis output directory
        """
        self.static_files = {key: value for key, value in self.expected_static_files().items()}

    def _init_files(self):
        """
        Initializes the paths to files that might not have consistent naming

        including: the targets .bed file with the chromosome target regions
        """
        self.set_file(name = 'targets_bed', path = find.find(search_dir = self.dir, inclusion_patterns = "*.bed", exclusion_patterns = '*.pad10.bed', search_type = 'file', num_limit = 1, level_limit = 0))

    def get_analysis_config(self):
        """
        Creates a dictionary of config values to pass to child Sample objects

        Returns
        -------
        dict
            a dictionary of configuration values
        """
        analysis_config = {}
        analysis_config['analysis_id'] = self.id
        analysis_config['analysis_dir'] = self.dir
        analysis_config['results_id'] = self.results_id

        analysis_config['dirs'] = self.dirs
        analysis_config['files'] = self.files
        analysis_config['static_files'] = self.static_files

        analysis_config['analysis_is_valid'] = self.is_valid

        # analysis_config['sns_config'] = self.sns_config
        return(analysis_config)

    def expected_static_files(self):
        """
        Creates a dictionary of files that are expected to exist in the analysis output
        Returns
        -------
        dict
            a dictionary of files that are expected to exist in the analysis dir
        """
        expected_files = {}
        # samplesheet file with the run's paired samples
        expected_files['paired_samples'] = os.path.join(self.dir, 'samples.pairs.csv')
        # file with the original starting .fastq file paths & id's
        expected_files['samples_fastq_raw'] = os.path.join(self.dir, 'samples.fastq-raw.csv')
        # file with settings for the analysis
        expected_files['settings'] = os.path.join(self.dir, 'settings.txt')
        # summary table produced at the end of the WES pipeline
        expected_files['summary_combined_wes'] = os.path.join(self.dir, 'summary-combined.wes.csv')
        return(expected_files)

    def get_qsub_logfiles(self, logdir = None):
        """
        Gets the list of log files from the analysis' qsub logs directory

        Parameters
        ----------
        logdir: str
            the path to the qsub log directory. If ``None``, a directory called ``logs-qsub`` will be searched for in the analysis output directory and used instead

        Returns
        -------
        list
            a list of file paths to qsub logs

        """
        log_files = []
        # try to get the logdir from self
        if not logdir:
            logdir = self.list_none(self.get_dirs('logs-qsub'))
        if not logdir:
            raise AnalysisItemMissing(message = 'Qsub log dir not found for the analysis', errors = '')
        else:
            # find all the log files
            for item in find.find(logdir, search_type = 'file'):
                log_files.append(item)
        return(log_files)

    def check_qsub_log_errors_present(self, log_files = None, err_patterns = ("ERROR:",)):
        """
        Checks the qsub log files for errors, by searching for lines that include known 'error' patterns

        Parameters
        ----------
        log_files: list
            a list of paths to qsub log files

        Returns
        -------
        bool
            ``True`` or ``False`` whether or not errors are detected in the qsub log files. If ``True``, the paths to files with error messages will be printed in the log.
        """
        contains_errors = {}
        # try to find the log files from self
        if not log_files:
            log_files = self.get_qsub_logfiles()
        if not log_files:
            raise AnalysisItemMissing(message = 'Qsub log files not found for the analysis.', errors = '')

        # check all the files for the patterns
        for log_file in log_files:
            with open(log_file, 'rb') as f:
                lines = f.readlines()
            for line in lines:
                for err_pattern in err_patterns:
                    if err_pattern in line:
                        contains_errors[log_file] = True

        # return a boolean for presence of errors
        if len(contains_errors) < 1:
            return(False)
        else:
            # True or False; any values are True = some log(s) contained error(s)
            if any(contains_errors.values()):
                self.logger.error('Error messages were found in some qsub logs')
                self.logger.debug('qsub log files containing errors: {0}'.format('\n'.join([path for path, value in contains_errors.items() if value == True])))
            return(any(contains_errors.values()))

    def get_summary_combined_contents(self, summary_combined_wes_file = None):
        """
        Gets the contents of the 'summary-combined.wes.csv' file

        Parameters
        ----------
        summary_combined_wes_file: str
            the path to the 'summary-combined.wes.csv' file, otherwise if ``None`` the file will be retrieved automatically from the analysis output

        Returns
        -------
        list
            a list of dictionaries representing the entries in the file, as read in by ``csv.DictReader``
        """
        # try to get the file from self if not passed
        if not summary_combined_wes_file:
            summary_combined_wes_file = self.static_files.get('summary_combined_wes', None)
        if not summary_combined_wes_file:
            raise AnalysisItemMissing(message = 'Could not find the summary_combined_wes_file ("summary-combined.wes.csv") for the analysis.', errors = '')

        # try to open it anyway
        rows = []
        with open(summary_combined_wes_file, 'rb') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                rows.append(row)
        return(rows)

    def summary_combined_contains_errors(self, summary_combined_wes_rows = None, err_pattern = 'X'):
        """
        Checks the 'summary-combined.wes.csv' file for errors; any entry in the sheet that looks like 'X'

        Parameters
        ----------
        summary_combined_wes_rows: list
            a list of dictionaries representing the entries in the file, as read in by ``csv.DictReader``
        err_pattern: str
            error pattern to search for in the file. Defaults to 'X'

        Returns
        -------
        bool
            ``True`` or ``False`` whether or not errors were detected in the file
        """
        # try to get the contents from self if not passed
        if not summary_combined_wes_rows:
            summary_combined_wes_rows = self.get_summary_combined_contents()
        if not summary_combined_wes_rows:
            # TODO: put an exception here
            self.logger.error('Could not find the summary_combined_wes_file contents')

        contains_errors = {}
        # try to parse it anyway
        for row in summary_combined_wes_rows:
            sampleID = row['#SAMPLE']
            # check all entries except the 'Sample'
            for item in [value for key, value in row.items() if key != '#SAMPLE']:
                if item == err_pattern:
                    contains_errors[sampleID] = True

        # return a boolean for presence of errors
        if len(contains_errors) < 1:
            return(False)
        else:
            # True or False; any values are True = some log(s) contained error(s)
            if any(contains_errors.values()): self.logger.warning('Error messages were found in "summary-combined.wes.csv" file for samples: {0}'.format([sampleID for sampleID, value in contains_errors.items() if value == True]))
            return(any(contains_errors.values()))

    def validate(self):
        """
        Checks if the analysis is considered valid for downstream usage

        Returns
        -------
        bool
            ``True`` or ``False`` whether or not the analysis passes validation criteria
        """
        self.validations = {}

        # make sure dir exists
        dir_validation = os.path.exists(self.dir)
        validation = {
            'dir_exists': {
            'status': os.path.exists(self.dir),
            'note': 'Whether or not the analysis directory ({0}) exists'.format(self.dir)
            }
        }
        self.validations.update(validation)


        # make sure all expected files exist
        expected_static_files_existences = [(key, value, os.path.exists(value)) for key, value in self.expected_static_files().items()]
        static_files_validations = {}
        validation = {
            'expected_static_files_exist': {
            'status': all([item[2] for item in expected_static_files_existences]),
            'note': 'Whether or not all of the expected files in the analysis exist;\n{0}'.format('\n'.join([str(i) for i in expected_static_files_existences]))
            }
        }
        self.validations.update(validation)

        # check for qsub log errors
        validation = {
            'no_qsub_log_errors_present': {
            'status': not self.check_qsub_log_errors_present(),
            'note': 'Whether or not errors are present in the qsub logs'
            }
        }
        self.validations.update(validation)

        # check for 'X' error entries in
        validation = {
            'no_summary_combined_errors': {
            'status': not self.summary_combined_contains_errors(),
            'note': 'Whether or not entries are present in the summary combined file'
            }
        }
        self.validations.update(validation)
        all_valid = [subdict['status'] for key, subdict in self.validations.items()]

        self.logger.debug('analysis validations:\n{0}'.format(json.dumps(self.validations, indent = 4)))

        is_valid = all(all_valid)
        self.logger.info('Analysis output passed validation: {0}'.format(is_valid))

        return(is_valid)


    def get_samplesIDs_from_samples_fastq_raw(self, samples_fastq_raw_file = None):
        """
        Gets the samples in the run from the samples_fastq_raw file

        Parameters
        ----------
        samples_fastq_raw_file: str
            path to the samples_fastq_raw_file for the analysis. If ``None``, it is retrieved automatically from the analysis output

        Returns
        -------
        list
            a list of sample ID's of the samples in the analysis
        """
        # self.logger.debug("Getting sample ID's from the 'samples_fastq_raw' file for the analysis")
        samplesIDs = []
        # try to get the file if it wasn't passed
        if not samples_fastq_raw_file:
            samples_fastq_raw_file = self.static_files.get('samples_fastq_raw', None)
        if samples_fastq_raw_file:
            with open(samples_fastq_raw_file, "rb") as csvfile:
                reader = csv.reader(csvfile)
                for row in reader:
                    samplesIDs.append(row[0])
        else:
            raise AnalysisItemMissing(message = 'The "samples_fastq_raw" file could not be found for the analysis.', errors = '')
        # unique entries only
        samplesIDs = list(set(samplesIDs))
        return(samplesIDs)

    def get_samples(self, samplesIDs = None):
        """
        Gets the samples for the analysis, as ``SnsAnalysisSample`` objects

        Parameters
        ----------
        samplesIDs: list
            a list of character strings representing sample ID's

        Returns
        -------
        list
            a list of ``SnsAnalysisSample`` objects
        """
        samples = []
        # try to get the sample IDs
        if not samplesIDs:
            samplesIDs = self.get_samplesIDs_from_samples_fastq_raw()
        for samplesID in samplesIDs:
            samples.append(SnsAnalysisSample(id = samplesID, analysis_config = self.get_analysis_config(), sns_config = self.sns_config, extra_handlers = self.extra_handlers))
        return(samples)





class SnsAnalysisSample(AnalysisItem):
    """
    Container for metadata about a single sample in the sns WES targeted exome sequencing run analysis output

    Examples
    --------
    Example usage::

        from sns_classes.classes import SnsWESAnalysisOutput
        import config
        config.sns.update({'email_recipients': 'foo@bar.edu'})
        d = '/ifs/data/molecpathlab/scripts/snsxt/snsxt/fixtures/sns_output/sns_analysis1'
        x = SnsWESAnalysisOutput(dir = d, id = 'sns_analysis1', sns_config = config.sns)
        samples = x.get_samples()
        sample = samples[0]
        sample.sns_config['analysis_output_index'].items()
        # note: this file retrieval method is deprecated and should not be used
        pattern = sample.id + '.dd.ra.rc.bam'
        sample.get_output_files(analysis_step = 'BAM-GATK-RA-RC', pattern = pattern)

    """

    def __init__(self, id, analysis_config, sns_config, extra_handlers = None):
        """
        Parameters
        ----------
        id: str
            ID for the sample
        analysis_config: dict
            dictionary of configuration settings passed on from the parent ``SnsWESAnalysisOutput`` object
        sns_config: dict
            dictionary of configuration settings passed on from the parent ``SnsWESAnalysisOutput`` object
        extra_handlers: list
            a list of Filehandlers, or ``None``
        """
        AnalysisItem.__init__(self, id = id, extra_handlers = extra_handlers)
        self.id = str(id)
        # set up per-sample logger
        # self.logger = log.build_logger(name = self.id)
        # if extra_handlers:
        #     self.logger = log.add_handlers(logger = self.logger, handlers = extra_handlers)
        # self.logger.debug("Initialized logging for sample: {0}".format(self.id))

        self.analysis_config = analysis_config
        self.sns_config = sns_config
        # self.logger.debug("Analysis is: {0}".format(self.analysis))

        # file matching pattern based on the sample's id
        self.search_pattern = '{0}*'.format(self.id)

        self._init_analysis_attrs()
        self._init_dirs()
        self._init_files()

    def _init_analysis_attrs(self, analysis_config = None):
        """
        Initializes object attributes passed from the parent analysis, for convenience

        Parameters
        ----------
        analysis_config: dict
            dictionary of configuration settings
        """
        if not analysis_config:
            analysis_config = self.analysis_config
        self.analysis_id = analysis_config['analysis_id']
        self.analysis_dir = analysis_config['analysis_dir']
        self.results_id = analysis_config['results_id']
        self.static_files = analysis_config['static_files']
        self.analysis_is_valid = analysis_config['analysis_is_valid']

    def _init_dirs(self, analysis_config = None):
        """
        Initializes the paths to dirs for the sample in the analysis

        Parameters
        ----------
        analysis_config: dict
            dictionary of configuration settings
        """
        if not analysis_config:
            analysis_config = self.analysis_config
        for name, paths in analysis_config['dirs'].items():
            if name not in ['_parent']:
                self.set_dir(name = name, path = paths)

    def _init_files(self, analysis_config = None):
        """
        Initializes the paths to files which were created dynamically based on the settings and samples of the analysis at run time

        Parameters
        ----------
        analysis_config: dict
            dictionary of configuration settings
        """
        if not analysis_config:
            analysis_config = self.analysis_config

        for name, paths in analysis_config['files'].items():
            self.set_file(name = name, path = paths)

    def get_output_files(self, analysis_step, pattern):
        """
        Gets a file from the sample's analysis output, based on the ``analysis_output_index`` config listing the expected file types at each output, in addition to the criteria specified by the function args

        Parameters
        ----------
        analysis_step: str
            the name of a directory in the analysis output from which to search for a sample's output file
        pattern: str
            a filename pattern to use when searching for the file

        Returns
        -------
        list
            a list of files for the sample from an analysis step
        """
        # get the dirpath for the analysis step from the analysis dir; return None if there isn't one set for the provided step
        search_dir = self.list_none(self.analysis_config['dirs'][analysis_step])
        patterns = [pattern, self.search_pattern]
        f = []
        if search_dir:
            # self.logger.debug("Searching for {0} files in {1}, dir: {2}".format(patterns, analysis_step, search_dir))
            f = find.find(search_dir = search_dir, inclusion_patterns = patterns, search_type = 'file', match_mode = 'all')
            # self.logger.debug('Found: {0}'.format(f))
        else:
            raise AnalysisItemMissing(message = "search_dir not found for {0}, dir: {1}".format(analysis_step, search_dir), errors = '')
        return(f)
