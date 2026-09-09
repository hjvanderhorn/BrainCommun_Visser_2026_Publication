import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

from statsmodels.stats.outliers_influence import variance_inflation_factor
from sklearn.metrics import roc_auc_score, log_loss
import statsmodels.formula.api as smf
import statsmodels.api as sm
from patsy import dmatrix

import warnings

# Treat warnings as exceptions
warnings.simplefilter("error", category=UserWarning)


def clean_concentrations(dataframe, biomarkers, biomarker_data):
    """
    codes unavailable data as 0; data below llod as 1/2 * llod; also below lloq as 1/2 lloq

    :param dataframe: pandas.dataframe, The dataframe holding the data
    :param biomarkers: list, An array of biomarkers which need to be processed
    :param biomarker_data: dict, {biomarker:llod}

    :return: None, manipulation occurs on array not a copy
    """

    # code unavailable V1 data as np.nan()
    dataframe.loc[dataframe['il6_conc'].isna(), ['il6_conc', 'il8_conc', 'il10_conc',
                                                 'nfl_conc', 'il6_conc_lloq', 'il8_conc_lloq', 'il10_conc_lloq']] = np.nan
    # dataframe.loc[dataframe['gfap_conc'].isna(), ['gfap_conc']] = np.nan
    dataframe.loc[dataframe['il6_conc_v2'].isna(), ['il6_conc_v2', 'il8_conc_v2', 'il10_conc_v2',
                                                    'nfl_conc_v2', 'il6_conc_v2_lloq', 'il8_conc_v2_lloq', 'il10_conc_v2_lloq']] = np.nan

    # Loop over the different biomarkers to set values LLOD as 1/2 LLOD
    for biomarker in biomarkers:

        dataframe[f'{biomarker}_v1_clean'] = dataframe[f'{biomarker}_conc'].apply(
            lambda x: biomarker_data[biomarker][1] if x == biomarker_data[biomarker][0] else x)
        dataframe[f'{biomarker}_v1_clean'] = dataframe[f'{biomarker}_v1_clean'].apply(
            lambda x: biomarker_data[biomarker][3] if x == biomarker_data[biomarker][2] else x)
        dataframe[f'{biomarker}_v1_clean'] = dataframe[f'{biomarker}_v1_clean'].astype(
            float)

        dataframe[f'{biomarker}_v2_clean'] = dataframe[f'{biomarker}_conc_v2'].apply(
            lambda x: biomarker_data[biomarker][1] if x == biomarker_data[biomarker][0] else x)
        dataframe[f'{biomarker}_v2_clean'] = dataframe[f'{biomarker}_v2_clean'].apply(
            lambda x: biomarker_data[biomarker][3] if x == biomarker_data[biomarker][2] else x)
        dataframe[f'{biomarker}_v2_clean'] = dataframe[f'{biomarker}_v2_clean'].astype(
            float)

        if biomarker in ['il6', 'il8', 'il10']:
            dataframe[f'{biomarker}_v1_lloq_clean'] = dataframe[f'{biomarker}_conc_lloq'].apply(
                lambda x: biomarker_data[biomarker][5] if x == biomarker_data[biomarker][4] else x)
            dataframe[f'{biomarker}_v1_lloq_clean'] = dataframe[f'{biomarker}_v1_lloq_clean'].astype(
                float)

            dataframe[f'{biomarker}_v2_lloq_clean'] = dataframe[f'{biomarker}_conc_v2_lloq'].apply(
                lambda x: biomarker_data[biomarker][5] if x == biomarker_data[biomarker][4] else x)
            dataframe[f'{biomarker}_v2_lloq_clean'] = dataframe[f'{biomarker}_v2_lloq_clean'].astype(
                float)
        else:
            continue

    return dataframe


def clean_poct_data(dataframe):
    dataframe['gfap_v1_clean'] = dataframe['gfap'].apply(
        lambda x: 30 if x == '<30' else x)
    dataframe['uchl1_v1_clean'] = dataframe['uchl1'].apply(
        lambda x: 200 if x == '<200' else x)
    dataframe['uchl1_v1_clean'] = dataframe['uchl1_v1_clean'].apply(
        lambda x: 3200 if x == '>3200' else x)
    dataframe['uchl1_v1_clean'].astype(float)


def createTimeGroup(minutes):
    """
    Groups time after trauma into 6 groups. 

    :param minutes: float | int, time in minutes since trauma

    :return: group, group which the given time in minutes belongs to.
    """

    if minutes <= 60:
        group = 1
    elif minutes <= 120:
        group = 2
    elif minutes <= 180:
        group = 3
    elif minutes <= 240:
        group = 4
    elif minutes <= 300:
        group = 5
    elif minutes > 300:
        group = 6
    else:
        return

    return group


def dichotimizeOutcome(gose):
    """
    Groups GOSE scores into complete recovery (gose = 8) or incomplete recovery (Gose < 8).

    :param gose: int | float, gose score six months after trauma

    :return: outcome, 0 = incomplete recovery, 1 = complete recovery, 99 = outcome unknown

    """

    if gose == 8:
        return 1
    elif gose < 8 and gose > 0:
        return 0
    else:
        return 99


def chooseTitle(biomarker):
    """ 
    Returns the full biomarker title from the abbreviation.

    :param biomarker: string, the abbreviated biomarker title

    :return: string, the unabbreviated biomarker title

    """

    if biomarker == 'il6':
        return 'Interleukin 6'
    elif biomarker == 'il8':
        return "Interleukin 8"
    elif biomarker == 'il10':
        return "Interleukin 10"
    elif biomarker == 'gfap':
        return "GFAP"
    elif biomarker == 'nfl':
        return "Neurofilament Light"
    elif biomarker == 'ft':
        return 'Free Thiols'


def createBoxPlot(biomarker, dataframe):
    """
    Create a boxplot for a particular biomarker of concentrations in the acute phase and subacute phase.

    :param biomarker: string, the biomarker to be modelled
    :param dataframe: pandas.dataframe, the dateframe that contains the data

    :return: matplotlib boxplot
    """
    fig, ax = plt.subplots(1, 2, sharey=True, dpi=300)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].yaxis.set_ticks_position('none')
    graph = sns.boxplot(data=dataframe.loc[dataframe[f'{biomarker}_v1_clean'] != 0], x='group', y=f'{biomarker}_v1_clean', order=[1, 0], ax=ax[0], boxprops=dict(
        facecolor='none', edgecolor='black', linewidth=.8), whiskerprops=dict(linewidth=.8), medianprops=dict(linewidth=.8), capprops=dict(linewidth=.8), showfliers=False, flierprops=dict(markeredgecolor='black'))
    ylimits = ax[0].get_ylim()
    sns.stripplot(data=dataframe.loc[dataframe[f'{biomarker}_v1_clean'] != 0], x='group', y=f'{biomarker}_v1_clean', order=[
                  1, 0], ax=ax[0], alpha=1, dodge='true', s=6, marker="$\circ$", jitter=0.3)
    ax[0].set_xlabel("Acute phase")
    graph = sns.boxplot(data=dataframe.loc[dataframe[f'{biomarker}_v2_clean'] != 0], x=dataframe.loc[dataframe['group'] != 2, 'group'], y=f'{biomarker}_v2_clean', order=[1, 0], ax=ax[1], boxprops=dict(
        facecolor='none', edgecolor='black', linewidth=.8), whiskerprops=dict(linewidth=.8), medianprops=dict(linewidth=.8), capprops=dict(linewidth=.8), showfliers=False, flierprops=dict(markeredgecolor='black'))
    sns.stripplot(data=dataframe.loc[dataframe[f'{biomarker}_v2_clean'] != 0], x='group', y=f'{biomarker}_v2_clean', order=[
                  1, 0], ax=ax[1], alpha=1, dodge='true', s=6, marker="$\circ$", jitter=0.3)

    print(ylimits)
    # ax[0].set(ylim=ylimits)
    ax[0].set(ylim=(-5.93, 3000))
    ax[1].set_xlabel("Subacute phase")
    ax[0].set(ylabel=f'[{chooseTitle(biomarker)}] pg/ml')
    ax[1].set(ylabel='')
    ax[0].grid(visible=False)
    ax[1].grid(visible=False)
    ax[0].set_xticklabels(['mTBI', 'HC1'])
    ax[1].set_xticklabels(['mTBI', 'HC1'])
    ax[0].xaxis.labelpad = 20
    ax[1].xaxis.labelpad = 20
    ax[0].legend([], [], frameon=False)
    ax[0].spines['left'].set_color('black')
    ax[0].spines['bottom'].set_color('black')
    ax[1].spines['bottom'].set_color('black')

    fig.tight_layout()
    # fig.savefig(f'{biomarker}.png', dpi=300, format='png', bbox_inches='tight')

    return [fig, ax]


def createTimePlot(dataframe, biomarker, tick_labels):
    """
    Given a biomarker plot its average concentration as a function of time since injury

    :param biomarker: string, The biomarker being modelled

    :return: matplotlib graph
    """

    fig, ax = plt.subplots(figsize=(10, 4), dpi=300)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sns.barplot(data=dataframe,
                x='time_to_v1',
                y=f'{biomarker}_v1_clean',
                n_boot=100,
                edgecolor=".5",
                linewidth=2,
                facecolor=(0, 0, 0, 0)
                )
    ax.set(ylabel=f'[{chooseTitle(biomarker)}] pg/ml')
    ax.set(xlabel='Time to sample, hours \n (n)')
    ax.grid(visible=False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    tick_locations = np.arange(6)
    new_labels = ["".join(x)
                  for x in zip(tick_labels[0::2], tick_labels[1::2])]
    plt.xticks(tick_locations, new_labels)
    ax.xaxis.labelpad = 20
    # ax.set_xticklabels(['0 - 1', '1 - 2', '2 - 3', '3 - 4', '4 - 5', '5 - 24'])

# Checking distributions of biomarker data and the presence of outliers.


def normalityCheck(dataframe_v1, dataframe_v2, biomarker, plot=False):

    print(f"Shapiro test for {biomarker} at timepoint 1:")
    print(stats.shapiro(dataframe_v1))
    print('')
    print(f"Shapiro test for {biomarker} at timepoint 2")
    print(stats.shapiro(dataframe_v2))

    if plot:

        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.hist(dataframe_v1, bins=20, edgecolor='k')
        plt.title(f'Histogram of {biomarker} at Timepoint 1')
        plt.xlabel('Value')
        plt.ylabel('Frequency')

        # Q-Q plot for timepoint 1
        plt.subplot(1, 2, 2)
        stats.probplot(dataframe_v1, dist="norm", plot=plt)
        plt.title(f'Q-Q Plot of {biomarker} at Timepoint 1')

        plt.show()

        #    Histogram for timepoint 2
        plt.figure(figsize=(12, 6))

        plt.subplot(1, 2, 1)
        plt.hist(dataframe_v2, bins=20, edgecolor='k')
        plt.title(f'Histogram of {biomarker} at Timepoint 2')
        plt.xlabel('Value')
        plt.ylabel('Frequency')

        # Q-Q plot for timepoint 2
        plt.subplot(1, 2, 2)
        stats.probplot(dataframe_v2, dist="norm", plot=plt)
        plt.title(f'Q-Q Plot of {biomarker} at Timepoint 2')

        plt.show()

    else:
        return


def normalityCheckSingle(dataframe_v1, biomarker):

    print(f"Shapiro test for {biomarker} at timepoint 1:")
    print(stats.shapiro(dataframe_v1))

    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.hist(dataframe_v1, bins=20, edgecolor='k')
    plt.title(f'Histogram of {biomarker} at Timepoint 1')
    plt.xlabel('Value')
    plt.ylabel('Frequency')

    # Q-Q plot for timepoint 1
    plt.subplot(1, 2, 2)
    stats.probplot(dataframe_v1, dist="norm", plot=plt)
    plt.title(f'Q-Q Plot of {biomarker} at Timepoint 1')

    plt.show()


def naturalLogTransform(dataframe, biomarkerArray):
    """
    Natural log transform a given biomarker concentration:

    :param dataframe: pandas.dataframe, the frame which contains the to use data
    :param biomarkerArray: list, list of all biomarkers which need to be log transformed
    :param biomarkerMetaData: list, list of dictionaries containing details about the biomarkers (LLOQ, LLOD)

    :return: new columns with the natural log transformed biomarker concentrations

    """

    for biomarker in biomarkerArray:

        if biomarker in ['gfap', 'uchl1']:
            dataframe[f'{biomarker}_v1_clean_ln'] = dataframe[
                f'{biomarker}_v1_clean'].apply(lambda x: np.cbrt(x))
            # dataframe[f'{biomarker}_v1_clean_ln_no'] = dataframe[f'{biomarker}_v1_clean_no'].apply(lambda x: np.log(x))
        else:
            dataframe[f'{biomarker}_v1_clean_ln'] = dataframe[
                f'{biomarker}_v1_clean'].apply(lambda x: np.log(x))
            dataframe[f'{biomarker}_v2_clean_ln'] = dataframe[
                f'{biomarker}_v2_clean'].apply(lambda x: np.log(x))
            # dataframe[f'{biomarker}_v1_clean_ln_no'] = dataframe[f'{biomarker}_v1_clean_no'].apply(lambda x: np.log(x))
            # dataframe[f'{biomarker}_v2_clean_ln_no'] = dataframe[f'{biomarker}_v2_clean_no'].apply(lambda x: np.log(x))

        # if biomarker in ['il6', 'il8', 'il10']:
        #     dataframe[f'{biomarker}_v1_lloq_clean_ln'] = dataframe[f'{biomarker}_v1_lloq_clean'].apply(lambda x: np.log(x) if x > biomarkerMetaData[biomarker][5] else x)
        #     dataframe[f'{biomarker}_v2_lloq_clean_ln'] = dataframe[f'{biomarker}_v2_lloq_clean'].apply(lambda x: np.log(x) if x > biomarkerMetaData[biomarker][5] else x)

    return


def logTransform(dataframe, biomarkerArray, suffix=['_v1_clean', '_v2_clean']):
    """
    Log transform a given biomarker concentration:

    :param dataframe: pandas.dataframe, the frame which contains the to use data
    :param biomarkerArray: list, list of all biomarkers which need to be log transformed
    :param biomarkerMetaData: list, list of dictionaries containing details about the biomarkers (LLOQ, LLOD)

    :return: new columns with the natural log transformed biomarker concentrations

    """

    for biomarker in biomarkerArray:

        if biomarker in ['gfap', 'uchl1']:
            dataframe[f'{biomarker}{suffix[0]}_log'] = dataframe[
                f'{biomarker}{suffix[0]}'].apply(lambda x: np.cbrt(x))
            # dataframe[f'{biomarker}_v1_clean_ln_no'] = dataframe[f'{biomarker}_v1_clean_no'].apply(lambda x: np.log(x))
        else:
            dataframe[f'{biomarker}{suffix[0]}_log'] = dataframe[
                f'{biomarker}{suffix[0]}'].apply(lambda x: np.log(x))
            dataframe[f'{biomarker}{suffix[1]}_log'] = dataframe[
                f'{biomarker}{suffix[1]}'].apply(lambda x: np.log(x))
            # dataframe[f'{biomarker}_v1_clean_ln_no'] = dataframe[f'{biomarker}_v1_clean_no'].apply(lambda x: np.log(x))
            # dataframe[f'{biomarker}_v2_clean_ln_no'] = dataframe[f'{biomarker}_v2_clean_no'].apply(lambda x: np.log(x))

        # if biomarker in ['il6', 'il8', 'il10']:
        #     dataframe[f'{biomarker}_v1_lloq_clean_ln'] = dataframe[f'{biomarker}_v1_lloq_clean'].apply(lambda x: np.log(x) if x > biomarkerMetaData[biomarker][5] else x)
        #     dataframe[f'{biomarker}_v2_lloq_clean_ln'] = dataframe[f'{biomarker}_v2_lloq_clean'].apply(lambda x: np.log(x) if x > biomarkerMetaData[biomarker][5] else x)

    return


def cuberootTransform(dataframe, biomarkerArray, suffix=['_v1_clean', '_v2_clean']):
    """
    Log transform a given biomarker concentration:

    :param dataframe: pandas.dataframe, the frame which contains the to use data
    :param biomarkerArray: list, list of all biomarkers which need to be log transformed
    :param biomarkerMetaData: list, list of dictionaries containing details about the biomarkers (LLOQ, LLOD)

    :return: new columns with the natural log transformed biomarker concentrations

    """

    for biomarker in biomarkerArray:

        if biomarker in ['gfap', 'uchl1']:
            dataframe[f'{biomarker}{suffix[0]}_log'] = dataframe[
                f'{biomarker}{suffix[0]}'].apply(lambda x: np.cbrt(x))
            # dataframe[f'{biomarker}_v1_clean_ln_no'] = dataframe[f'{biomarker}_v1_clean_no'].apply(lambda x: np.log(x))
        else:
            dataframe[f'{biomarker}{suffix[0]}_cbrt'] = dataframe[
                f'{biomarker}{suffix[0]}'].apply(lambda x: np.cbrt(x))
            dataframe[f'{biomarker}{suffix[1]}_cbrt'] = dataframe[
                f'{biomarker}{suffix[1]}'].apply(lambda x: np.cbrt(x))
            # dataframe[f'{biomarker}_v1_clean_ln_no'] = dataframe[f'{biomarker}_v1_clean_no'].apply(lambda x: np.log(x))
            # dataframe[f'{biomarker}_v2_clean_ln_no'] = dataframe[f'{biomarker}_v2_clean_no'].apply(lambda x: np.log(x))

        # if biomarker in ['il6', 'il8', 'il10']:
        #     dataframe[f'{biomarker}_v1_lloq_clean_ln'] = dataframe[f'{biomarker}_v1_lloq_clean'].apply(lambda x: np.log(x) if x > biomarkerMetaData[biomarker][5] else x)
        #     dataframe[f'{biomarker}_v2_lloq_clean_ln'] = dataframe[f'{biomarker}_v2_lloq_clean'].apply(lambda x: np.log(x) if x > biomarkerMetaData[biomarker][5] else x)

    return


def getColumn(n_variables, current_index, row):
    if row == 1:
        if current_index == (n_variables / 2):
            column = 0
        elif current_index == (n_variables / 2) + 1:
            column = 1
        elif current_index == (n_variables / 2) + 2:
            column = 2
    else:
        column = current_index
    return column


def getRow(case, current_index, n_variables):
    print(case)
    if (case == 1) & (current_index >= (n_variables/2)):
        row = 1
        column = getColumn(n_variables, current_index, row)
    elif (case == 1) & (current_index < (n_variables / 2)):
        row = 0
        column = getColumn(n_variables, current_index, row)
    elif case == 2:
        row = 0
        column = current_index
    elif (case == 3) & ((current_index >= (n_variables // 2) + 1)):
        row = 1
        column = getColumn(n_variables + 1, current_index, row)
    else:
        row = 0
        column = getColumn(n_variables, current_index, row)

    return row, column


def check_distribution(dataframe, variables):
    n_variables = len(variables)

    if n_variables > 1:
        if n_variables < 4:
            fig, ax = plt.subplots(1, n_variables, figsize=(10, 6), dpi=300)
            case = 2
        elif (n_variables % 2) == 0:
            fig, ax = plt.subplots(
                2, int(n_variables / 2), figsize=(10, 10), dpi=300)
            case = 1

        else:
            fig, ax = plt.subplots(
                2, int((n_variables // 2) + 1), figsize=(10, 10), dpi=300)
            case = 3

        for index, variable in enumerate(variables):
            row, column = getRow(case, index, n_variables)
            int_data = dataframe.loc[dataframe[variable].notna()]
            shapiro = stats.shapiro(int_data[variable])
            ks = stats.kstest(int_data[variable], stats.norm.cdf)
            print(variable, ':')
            print(f"Shapiro: {shapiro}")
            print(f"Kolmogorov_Smirnov: {ks}")
            if case == 2:
                sns.histplot(data=int_data, x=variable, ax=ax[column])
            else:
                sns.histplot(data=int_data, x=variable, ax=ax[row, column])

    else:
        int_data = dataframe.loc[dataframe[variables[0]].notna()]
        shapiro = stats.shapiro(int_data[variables[0]])
        ks = stats.kstest(rvs=int_data[variables], cdf=stats.norm.cdf)
        print(variables[0], ':')
        print(f"Shapiro: {shapiro}")
        print(f"Kolmogorov_Smirnov: {ks}")
        sns.histplot(data=int_data, x=variables[0])
        return


def getCorrelationValue(dataframe, parameter1, parameter2, correlation_stat):
    data_nona = dataframe.loc[(dataframe[parameter1].notna()) & (
        dataframe[parameter2].notna())]

    if correlation_stat == 'p':
        result = stats.pearsonr(data_nona[parameter1], data_nona[parameter2])
        print(f"Pearsons correlation is r"
              f"{round(result.statistic, 4)}, significance {round(result.pvalue, 4)}")
    elif correlation_stat == 's':
        result = stats.spearmanr(data_nona[parameter1], data_nona[parameter2])
        print(f"Spearmans rho: "
              f"{round(result.statistic, 4)}, significance {round(result.pvalue, 4)}")
    else:
        print(f"Provided correlation method doesn't exits {correlation_stat}")


def createZScores(value, std, mean):

    z_score = (value - mean) / std
    return z_score


def createColumnNoOutlier(z_criteria, metabolite, dataframe):
    metabolite_z = f"{metabolite}_z"

    std = np.std(dataframe[metabolite])
    mean = np.mean(dataframe[metabolite])

    dataframe[metabolite_z] = dataframe[metabolite].apply(
        lambda x: createZScores(x, std, mean))
    outliers = dataframe.loc[(dataframe[metabolite_z] >= z_criteria) | (
        dataframe[metabolite_z] <= -z_criteria), [metabolite, metabolite_z, 'record_id', 'group']]
    print(outliers)

    dataframe[f"{metabolite}_no"] = dataframe[[metabolite_z, metabolite]].apply(
        lambda x: x[1] if (x[0] < z_criteria) and (x[0] > -z_criteria) else np.nan, axis=1)


def calculate_extracranial(ais_head, iss):

    if (np.isnan(ais_head)) or (np.isnan(iss)):
        return 99

    extracranial = 0
    if (ais_head == 1) & (iss > 2):
        extracranial = 1
    elif (ais_head == 2) & (iss > 5):
        extracranial = 1
    elif (ais_head == 3) & (iss > 10):
        extracranial = 1
    elif (ais_head == 4) & (iss > 17):
        extracranial = 1
    elif (ais_head == 0) & (iss > 0):
        extracranial = 1

    return extracranial


def getSpecDataString(dataframe, variable, variable_data):

    dataStringArray = []

    if (variable == 'gcs') or (variable == 'centre') or (variable == 'occupation_status') or (variable == 'education_type'):
        for subvariable in variable_data[1]:
            dataString, _ = getCharacteristicString(
                dataframe, variable, ['ord', [subvariable]])
            dataStringArray.append(dataString)

    elif variable == 'injury_mechanism':
        injury_mechanism_groups = [[1, 3, 8, 9], [2], [4], [5], [6, 7]]
        for injury_group in injury_mechanism_groups:
            total = 0
            for mechanism in injury_group:
                dataString, _ = getCharacteristicString(
                    dataframe, variable, ['ord', [mechanism]])
                amount = int(dataString[0:2])
                total = total + amount

            denominator = len(dataframe[variable])
            percentageOfTotal = round((total / denominator) * 100, 2)
            newDataString = f"{total} ({percentageOfTotal}%)"
            dataStringArray.append(newDataString)

    return dataStringArray


def getCharacteristicString(dataframe, variable, variable_data):
    dataType = variable_data[0]

    total = len(dataframe.loc[dataframe[variable].notna()])

    if dataType == 'cont':
        mean = round(dataframe[variable].mean(), 0)
        std = round(dataframe[variable].std(), 0)
        min = round(dataframe[variable].min(), 0)
        max = round(dataframe[variable].max(), 0)

        dataString = f"{mean} {std} ({min}-{max})"

    if dataType == 'likert':
        median = dataframe[variable].median()
        lq = dataframe[variable].quantile(0.25)
        uq = dataframe[variable].quantile(0.75)

        dataString = f"{median} [{lq}, {uq}]"

    if dataType == 'ord':
        numerator = variable_data[1][0]
        numerator_data = dataframe.loc[dataframe[variable] == numerator, [
            variable]]

        count = len(numerator_data)
        denominator = len(dataframe[variable])
        percentage = round((count / denominator) * 100, 1)

        dataString = f"{count} ({percentage}%)"

    if dataType == 'spec':

        dataString = getSpecDataString(dataframe, variable, variable_data)

    return [dataString, total]


def createCharacteristicTable(variables, groups, groupDataframes):

    tableData = {'CR': [], 'IR': []}
    rowLabels = ['count']
    tableData['CR'].append(len(groupDataframes['CR']))
    tableData['IR'].append(len(groupDataframes['IR']))

    totals = {
        'CR': [],
        'IR': [],
        'variables': []
    }

    for variable_ in variables:
        totals['variables'].append(variable_)
        for group in groups:
            information, total = getCharacteristicString(
                groupDataframes[group], variable_, variables[variable_])
            totals[group].append(total)

            if variables[variable_][0] == 'spec':
                for i in range(len(variables[variable_][2])):
                    subvariable = variables[variable_][2][i]
                    if subvariable not in rowLabels:
                        rowLabels.append(subvariable)

                    tableData[group].append(information[i])
            else:
                if variable_ not in rowLabels:
                    rowLabels.append(variable_)
                tableData[group].append(information)

    table_size = len(tableData['CR'])

    # for i in range(2, table_size):
    #     tableData['hc2'][i] = '-'
    #     tableData['hc1'][i] = '-'

    # print(len(tableData['hc2']))
    tableData['Variables'] = rowLabels

    characteristicTable = pd.DataFrame.from_dict(tableData, orient='columns')
    characteristicTable.set_index('Variables', inplace=True)

    totals_dataframe = pd. DataFrame.from_dict(totals)
    totals_dataframe.set_index('variables', inplace=True)

    display(totals_dataframe)

    return characteristicTable
