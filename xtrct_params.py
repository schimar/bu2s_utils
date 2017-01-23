#! /usr/bin/python
#
# This script reads through a bu2s parameters file and prints the contents (multi-row) in one row. # used with a shell script (called lsBUnzip2paramS.sh), which loops through folders (of individual simulation runs), bunzips the paremeters.m.bz2 file, calls this python script, and bzips the file again. The header info can be found at the bottom of this script.
# NOTE: 1) the script will print a warning, without adding the current run's parameters to the output file, if any bu2s version other then 3.6.1 was used; 2) if a parameter can have multiple values, these are output as space-separated, whereas columns are tab-delimited

# Usage: ./xtrc_params.py parameters.m <RunFolderName>

from sys import argv

runName = argv[2]

with open(argv[1], 'r') as file:
    #hdrs = dict()
    newline = list()
    newline.append(runName)
    for lineORG in file:
        line = lineORG.split('\n')[0].split(' = ')
        if line[0] == 'CodeVersion':
            if line[1] != "'bu2s_3.6.1';":
                print 'Warning: %s is version %s; parameters not added' % (runName, line[1])
                #print type(line[1])
            else:
                continue
        elif line[0] == '':
            continue
        elif lineORG[0] == '%':
            continue
        else:
            if line[1][0] == '[':
                value = line[1].split('[')[1].split(']')[0]
            else:
                value = line[1].split(';')[0]
            newline.append(value)
            #print line[0]
    print '\t'.join(newline)

    file.close()





# headers (column names) for bu2s version 3.6.1

#D	PATCHES	nPATCHES	USE_MUTATIONS_FROM_FILES	ADDITIVE_FITNESS	MULTIPLICATIVE_FITNESS	TWO_DEME	DETERMINISTIC	nGENERATIONS_MAX	nMUTATIONS	INITIAL_CONDITIONS	INITIAL_POPULATION_SIZE	MEAN_S	SYMMETRIC_MUTATIONS	MAP_TYPE	GAMETE_PRODUCTION_MODE	FIXED_N	K	MOSAIC	SD_MOVE	OFFSPRING_IN_RANDOM_LOCATIONS	nCHROMOSOMES	TOTAL_MAP_LENGTH	MUTATION_DISTRIBUTION	START_WITH_BIG_IN_PLACE	BIG_START_S_VAL	FATTEN_TAIL_PROPORTION	FATTEN_TAIL_MAX	NUMBER_BIG_TO_ESTABLISH	END_PERIOD_ALLOPATRY	START_PERIOD_ALLOPATRY	PERIOD_ALLOPATRY	DMI_MEAN_EFFECT	DMI_PROBABILITY	PROBABILITY_DERIVED_REVERSED	POSITIVE_EPI_PROBABILITY	MEAN_POSITIVE_EPI_COEFF	START_THRESH_FOR_LD	END_THRESH_FOR_LD	MUTATIONS_PER_GENERATION	DEME0_CONSTANT_S	DEME1_CONSTANT_S	MU	FRACTION_SELECTED_LOCI	PURE_NEUTRAL_MODEL	PURE_ALLOPATRY_MODEL	epi_patch_multipliers	H	TS_SAMPLING_FREQUENCY	RECORDING_TIMES_IN_FILE	RI_THRESH	RI_REACHED	RECORD_LD_VALUES	BeginRecordingLD	LD_LOCI_SUBSAMPLE	LD_LowerBound	FST_MIN_RECORDING_THRESH	approachingSpeciationThreshold	firstTimeResNearMaxFit	lastTimeResNearImm	lastTimeResNearRand	totalGenerationsElapsed	estRunLength	totalMutationsIntroduced	N	nLOCI	nVariableLoci	totalDMIs	totalPositiveEpis	CONSIDER_EPISTASIS	LOCI_PER_CHROMOSOME
