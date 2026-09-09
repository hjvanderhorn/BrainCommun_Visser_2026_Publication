import warnings
import statsmodels.formula.api as smf
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from scikits.bootstrap import ci
from sklearn.metrics import roc_auc_score, log_loss
from tabulate import tabulate
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=RuntimeWarning,
                        message="overflow encountered in exp")


class logit_regression_model:
    """
    A logistic regression model wrapper for fitting, evaluating, and summarizing results.

    Attributes:
        predictors (list): List of predictor variable names.
        dv (str): Name of the dependent (response) variable.
        data (pd.DataFrame): Processed dataframe with predictors and dependent variable.
        X (pd.DataFrame): Predictor variables used for modeling.
        y (pd.Series): Response variable.
        predictor_string (str): Formula string for logistic regression.
        rsquare_comparison_model (float): Baseline R-squared value for comparison.
        c_stat_comparison_model (float): Baseline C-statistic value for comparison.
    """

    def __init__(self, predictors, dv, reg_data, interactions, r_square_base, c_stat_base):
        """
        Initializes the LogitRegressionModel with predictors, response variable, and comparison metrics.

        Args:
            predictors (list): Predictor variables.
            dv (str): Dependent variable.
            reg_data (pd.DataFrame): The input dataset.
            interactions (dict): Interaction terms to include in the model.
            r_square_comparison (float): R-squared value for comparison.
            c_stat_comparison (float): C-statistic value for comparison.
        """
        self.predictors = predictors
        self.dv = dv
        self.data = reg_data[predictors + [dv]].dropna(axis=0)
        self.X = self.data[predictors]
        self.y = self.data[dv]
        self.predictor_string = get_predictor_string(predictors, interactions)
        self.nagelkerke_base_model = r_square_base
        self.c_stat_base_model = c_stat_base

    def fit_main_model(self):
        formula = f"{self.dv} ~ {self.predictor_string}"
        print(formula)
        self.model = smf.logit(
            formula, data=self.data).fit(method="newton", disp=0, maxiter=1000)
        llf = self.model.llf
        self.deviance = -2 * llf
        self.dof = len(self.predictors)

    def get_model_summary(self):
        print(self.model.summary())

    def calculate_bootstrap_base(self, *bootstrap_indices, predictor_string):
        try:
            current_index_values = list(bootstrap_indices[0])
            new_data = get_dataframe_from_multiple_index(
                self.data, current_index_values, self.predictors, self.dv)
            # Scale the data

            boot_model = smf.logit(
                f"{self.dv} ~ {predictor_string}", data=new_data).fit(method="newton", disp=0, maxiter=1000)

            # Raw predicted values from original data and fitted model
            y_actual = new_data[self.dv]
            y_pred_raw = boot_model.predict()

            nagelkerke_rsquared = calculate_nagelkerke_new(
                y_actual, y_pred_raw)

            c_stat = roc_auc_score(y_actual, y_pred_raw)

            return nagelkerke_rsquared, c_stat

        except Exception as e:
            # print(f"bootstrap iteration failed: {e}")
            return np.nan, np.nan

    def calculate_bootstrap_delta(self, *bootstrap_indices, predictor_string_base, predictor_string_comp):
        try:
            current_index_values = list(bootstrap_indices[0])
            new_data = get_dataframe_from_multiple_index(
                self.data, current_index_values, self.predictors, self.dv)

            # Calculating statistics on the base model for this bootstrap iteration

            boot_model_base = smf.logit(
                f"{self.dv} ~ {predictor_string_base}", data=new_data).fit(method="newton", disp=0, maxiter=1000)

            y_actual = new_data[self.dv]
            y_pred_raw_base = boot_model_base.predict()

            nagelkerke_rsquared_base = calculate_nagelkerke_new(
                y_actual, y_pred_raw_base)

            c_stat_base = roc_auc_score(y_actual, y_pred_raw_base)

            # Calculating statistics on the comparison (delta) model for this bootstrap iteration

            boot_model_comp = smf.logit(
                f"{self.dv} ~ {predictor_string_comp}", data=new_data).fit(method="newton", disp=0, maxiter=1000)

            y_pred_raw_comp = boot_model_comp.predict()

            nagelkerke_rsquared_comp = calculate_nagelkerke_new(
                y_actual, y_pred_raw_comp)

            c_stat_comp = roc_auc_score(y_actual, y_pred_raw_comp)

            # calculating the delta statistics

            nagelkerke_rsquared_delta = nagelkerke_rsquared_comp - nagelkerke_rsquared_base
            c_stat_delta = c_stat_comp - c_stat_base

            return nagelkerke_rsquared_comp, c_stat_comp, nagelkerke_rsquared_delta, c_stat_delta
        except Exception as e:
            # print(f"bootstrap iteration failed: {e}")
            return np.nan, np.nan, np.nan, np.nan

    def calculate_bootstrap_optimism(self, *data, predictor_string_base, predictor_string_comp):
        try:
            current_index_values = list(data[0])
            new_data = get_dataframe_from_multiple_index(
                self.data, current_index_values, self.predictors, self.dv)

            y_actual = new_data[self.dv]
            # Calculating statistics on the base model for this bootstrap iteration

            boot_model_base = smf.logit(
                f"{self.dv} ~ {predictor_string_base}", data=new_data).fit(method="newton", disp=0, maxiter=1000)

            y_pred_raw_base = boot_model_base.predict()

            apparent_nagelkerke_rsquared_base = calculate_nagelkerke_new(
                y_actual, y_pred_raw_base)
            apparent_c_stat_base = roc_auc_score(y_actual, y_pred_raw_base)

            # Calculating statistics on the comparison (delta) model for this bootstrap iteration
            boot_model_comp = smf.logit(
                f"{self.dv} ~ {predictor_string_comp}", data=new_data).fit(method="newton", disp=0, maxiter=1000)

            y_pred_raw_comp = boot_model_comp.predict()
            apparent_nagelkerke_rsquared_comp = calculate_nagelkerke_new(
                y_actual, y_pred_raw_comp)
            apparent_c_stat_comp = roc_auc_score(y_actual, y_pred_raw_comp)

            # calculating the delta statistics
            apparent_nagelkerke_rsquared_delta = apparent_nagelkerke_rsquared_comp - \
                apparent_nagelkerke_rsquared_base
            apparent_c_stat_delta = apparent_c_stat_comp - apparent_c_stat_base

            # Using the newly created models to predict on the original data
            nagelkerke_r_squared_base_original = calculate_nagelkerke_new(
                self.y, boot_model_base.predict(self.X))
            nagelkerke_r_squared_comp_original = calculate_nagelkerke_new(
                self.y, boot_model_comp.predict(self.X))
            c_stat_base_original = roc_auc_score(
                self.y, boot_model_base.predict(self.X))
            c_stat_comp_original = roc_auc_score(
                self.y, boot_model_comp.predict(self.X))

            optimism_nagelkerke_original = apparent_nagelkerke_rsquared_base - \
                nagelkerke_r_squared_base_original
            optimism_nagelkerke_comp = apparent_nagelkerke_rsquared_comp - \
                nagelkerke_r_squared_comp_original

            delta_optimism_nagelkerke = optimism_nagelkerke_comp - optimism_nagelkerke_original

            optimism_c_stat_original = apparent_c_stat_base - c_stat_base_original
            optimism_c_stat_comp = apparent_c_stat_comp - c_stat_comp_original

            delta_optimism_c_stat = optimism_c_stat_comp - optimism_c_stat_original

            # saving optimsm corrected values for current iteration to ensure we get a distribution of optimsm corrected metrics
            corrected_nagelkerke_r2 = apparent_nagelkerke_rsquared_comp - optimism_nagelkerke_comp
            corrected_c_stat = apparent_c_stat_comp - optimism_c_stat_comp
            corrected_delta_r2 = apparent_nagelkerke_rsquared_delta - delta_optimism_nagelkerke
            corrected_delta_c = apparent_c_stat_delta - delta_optimism_c_stat

            return (apparent_nagelkerke_rsquared_comp,  # 1 nagelkerke of comparison model
                    apparent_c_stat_comp,  # 2 c-stat of comparison model
                    # 3 delta nagelkere of comparison model to reference model
                    apparent_nagelkerke_rsquared_delta,
                    apparent_c_stat_delta,  # 4 delta c-stat of comparison model to reference model
                    # 5 optimism of delta nagelkerke on bootstrap vs original data
                    delta_optimism_nagelkerke,
                    delta_optimism_c_stat,  # 6 optimism of delta c-stat on bootstrap vs original data
                    # 7 corrected nagelkerke of comparison model on current bootstrap iteration
                    corrected_nagelkerke_r2,
                    corrected_c_stat,  # 8 corrected c-stat of comparison model on current bootstrap iteration
                    # 9 corrected delta nagelkerke of comparison model to reference model on current bootstrap iteration
                    corrected_delta_r2,
                    # 10 corrected delta c-stat of comparison model to reference model on current bootstrap iteration
                    corrected_delta_c)

        except ValueError as e:
            # print(f"bootstrap iteration failed (value error): {e}")
            return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

        except Exception as e:
            # print(f"bootstrap iteration failed: {e}")
            return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    def calculate_nagelkerke(self):
        y_pred = self.model.predict()
        self.nagelkerke_rsquared = calculate_nagelkerke_new(self.y, y_pred)

    def calculate_c_stat(self):
        y_pred_raw = self.model.predict()
        logit_roc_auc = roc_auc_score(self.y, y_pred_raw)

        self.c_stat = logit_roc_auc

    def calculate_delta_nagelkerke(self):
        delta_nagelkerke = self.nagelkerke_rsquared - self.nagelkerke_base_model
        self.delta_nagelkerke = delta_nagelkerke

    def calculate_delta_c_stat(self):
        delta_c_stat = self.c_stat - self.c_stat_base_model
        self.delta_c_stat = delta_c_stat

    def calculate_all_stats(self):
        self.calculate_nagelkerke()
        self.calculate_c_stat()
        self.calculate_delta_nagelkerke()
        self.calculate_delta_c_stat()

    def bootstrap_statistic(self, statistic, n_iterations, ci_calculation_method="percentile", predictor_string=None):

        index_list = (get_index_list_dataframe(self.data), )

        predictor_string_base = None

        if predictor_string != None:
            predictor_string_base = predictor_string
            predictor_string_comp = self.predictor_string

        else:
            predictor_string_base = self.predictor_string

        if statistic == 'base_model':
            result = stats.bootstrap(
                data=index_list,
                statistic=lambda data: self.calculate_bootstrap_base(
                    data, predictor_string=predictor_string_base),
                n_resamples=n_iterations,
                axis=0,
                method=ci_calculation_method
            )

            self.ci_base = result

        elif statistic == 'delta':
            result = stats.bootstrap(
                data=index_list,
                statistic=lambda data: self.calculate_bootstrap_delta(
                    data, predictor_string_base=predictor_string_base, predictor_string_comp=predictor_string_comp),
                n_resamples=n_iterations,
                axis=0,
                method=ci_calculation_method
            )

            self.ci_delta = result

        elif statistic == 'delta_optimism':
            result = stats.bootstrap(
                data=index_list,
                statistic=lambda data: self.calculate_bootstrap_optimism(
                    data, predictor_string_base=predictor_string_base, predictor_string_comp=predictor_string_comp),
                n_resamples=n_iterations,
                axis=0,
                method=ci_calculation_method
            )

            self.ci_delta_optimism = result

# Utility functions wrapping the logreg model


def fitModel(predictors, dv, dataframe, isBaseModel, base_model_r2, base_model_c, summary=False, interactions={}):

    if isBaseModel:
        model = logit_regression_model(
            predictors, dv, dataframe, interactions, None, None)
        model.fit_main_model()
    else:
        # In this case it is a comparison model, therefore we provide r2 and c-stat form base model
        model = logit_regression_model(
            predictors, dv, dataframe, interactions, base_model_r2, base_model_c)
        model.fit_main_model()

    if summary:
        print('-'*80)
        if isBaseModel:
            (print('-> Baseline Model'))
        else:
            print('-> Comparison Model')
        model.get_model_summary()
        print('')

    return model


def calculate_statistics(model, display=True, CI=False, qc=False, n_iterations=1000, ci_calculation_method='percentile'):
    model.calculate_nagelkerke()
    model.calculate_c_stat()
    if CI:
        model.bootstrap_statistic(
            statistic='base_model', n_iterations=n_iterations, ci_calculation_method=ci_calculation_method)

    if qc:
        qc_bootstrap(model)

    if display:
        display_statistics(model, CI)


def display_statistics(model_class, CI=False):
    if CI:
        ci_nagelkerke = calculate_ci(
            model_class.ci_base.bootstrap_distribution[0])

        nagelkerke_bootstrap_ci = [
            ci_nagelkerke[0], ci_nagelkerke[1]]
        nagelkerke_score = model_class.nagelkerke_rsquared

        print(f"R2 = {round(nagelkerke_score, 4)}  "
              f"[{round(nagelkerke_bootstrap_ci[0], 4)} - {round(nagelkerke_bootstrap_ci[1], 4)}]")

        ci_c_stat = calculate_ci(model_class.ci_base.bootstrap_distribution[1])

        c_stat_bootstrap_ci = [
            ci_c_stat[0], ci_c_stat[1]]
        c_stat = model_class.c_stat

        print(f"C-statistic = {round(c_stat, 4)}  "
              f"[{round(c_stat_bootstrap_ci[0], 4)} - {round(c_stat_bootstrap_ci[1], 4)}]")

    else:
        nagelkerke_score = model_class.nagelkerke_rsquared
        print(f"Nagelkerke R2 = {round(nagelkerke_score, 4)}")

        c_stat = model_class.c_stat
        print(f"C-stat = {round(c_stat, 4)}")


def qc_bootstrap(model_class, type='base'):
    if type == 'base':
        bootstrap_array = model_class.ci_base.bootstrap_distribution
        columns = ['bootstrap_r2', 'bootstrap_c']
    elif type == 'delta':
        bootstrap_array = model_class.ci_delta.bootstrap_distribution
        columns = ['bootstrap_r2_delta', 'bootstrap_c_delta']
    else:
        bootstrap_array = model_class.ci_delta_optimism.bootstrap_distribution
        columns = ['bootstrap_r2_delta_op', 'bootstrap_c_delta_op']

    bootstrap_df = pd.DataFrame(bootstrap_array.T, columns=columns)

    create_hist_boot_distribution(bootstrap_array[0], columns[0])
    create_hist_boot_distribution(bootstrap_array[1], columns[1])

    return


# Next set of wrappers will be for comparison bootstraps

def compareModels(baseModel, predictors_comparison, dv, dataframe, spline_data, dataframe_spline, compare_models=True, verbose=False, interactions={}):
    if spline_data:
        compModel = fitModel(predictors_comparison, dv, dataframe_spline, False,
                             baseModel.nagelkerke_rsquared, baseModel.c_stat, verbose, interactions)
        datatype = 'spline_data'
    else:
        compModel = fitModel(
            predictors=predictors_comparison,
            dv=dv,
            dataframe=dataframe,
            isBaseModel=False,
            base_model_r2=baseModel.nagelkerke_rsquared,
            base_model_c=baseModel.c_stat,
            summary=True, interactions=interactions)
        datatype = 'models'

    # initializes nagelkerke and c-stat for the comparison model. Can be displayed as is or used for bootstrapping
    compModel.calculate_all_stats()

    if compare_models:
        # display_statistics_adjusted(compModel)

        print('')
        lr_pvalue = likelihood_ratio_test(
            baseModel.deviance, compModel.deviance, baseModel.dof, compModel.dof, datatype)

    return compModel, lr_pvalue


# Simple comparison wrapper i.e. no bootstrap necessary

def calculate_statistic_delta(baseModel, compModel, n_iterations, ci_calculation_method, biomarker, output_file_name):

    compModel.bootstrap_statistic(statistic="delta_optimism",
                                  n_iterations=n_iterations,
                                  predictor_string=baseModel.predictor_string,
                                  ci_calculation_method=ci_calculation_method)

    display_statistics_delta(compModel, biomarker, output_file_name)

    return


def simpleCompare(predictors_base, predictors_comp, dv, data, interactions, n_iterations=1000, ci_calculation_method='bca'):

    baseModel = fitModel(predictors_base, dv, data, True,
                         None, None, True)

    calculate_statistics(model=baseModel, display=True, CI=True)

    predictors_comp = predictors_comp

    compModel, lr_pvalue = compareModels(baseModel, baseModel.predictors + predictors_comp, dv, data, spline_data=False,
                                         dataframe_spline=None, compare_models=True, verbose=False, interactions=interactions)

    calculate_statistic_delta(baseModel, compModel,
                              n_iterations, ci_calculation_method, predictors_comp, 'test_1_24-06')

    return


def calculate_ci(boot):
    check_ci = np.nanpercentile(boot, [2.5, 97.5])
    return check_ci


def display_statistics_delta(compModel, biomarker, output_file_name, qc=False):

    output = {'R2_comp': [], 'c_stat_comp': [],
              'R2_delta': [], 'c_stat_delta': []}
    keys = ['R2_comp', 'c_stat_comp', 'R2_delta', 'c_stat_delta']
    numbers = [0, 1, 8, 9]

    count = 0

    for number in numbers:
        ci_check = calculate_ci(
            compModel.ci_delta_optimism.bootstrap_distribution[number])
        output[keys[count]].append(ci_check)
        count += 1

    nagelkerke_delta_optimism = np.nanmean(
        compModel.ci_delta_optimism.bootstrap_distribution[4])

    c_stat_delta_optimism = np.nanmean(
        compModel.ci_delta_optimism.bootstrap_distribution[5])

    if qc:
        create_hist_boot_distribution(
            compModel.ci_delta_optimism.bootstrap_distribution[4], 'nagelkerke_R2')

        create_hist_boot_distribution(
            compModel.ci_delta_optimism.bootstrap_distribution[5], 'c_stat')

    data = [
        ["R2", round(compModel.nagelkerke_rsquared, 4), np.nan, round(output['R2_comp'][0][0], 4), round(
            output['R2_comp'][0][1], 4)],
        ["C", round(compModel.c_stat, 4), np.nan, round(
            output['c_stat_comp'][0][0], 4), round(output['c_stat_comp'][0][1], 4)],
        ["Delta: R2", round(compModel.delta_nagelkerke, 4), nagelkerke_delta_optimism, round(
            output['R2_delta'][0][0], 4), round(output['R2_delta'][0][1], 4)],
        ["Delta: C", round(compModel.delta_c_stat, 4), c_stat_delta_optimism, round(
            output['c_stat_delta'][0][0], 4), round(output['c_stat_delta'][0][1], 4)],
    ]

    # Tabulate
    headers = ["Statistic", "Value", "Optimism", "CI Lower", "CI Upper"]
    table = tabulate(data, headers=headers, tablefmt="pretty")

    print(table)

    with open(f"./outputs/{output_file_name}.txt", 'a') as f:
        # Print the table
        print(f'Output summary for {biomarker}', file=f)
        print(table, file=f)

# Returns a list of all the indices of a dataframe. This can be used in bootstrapping methods


def get_index_list_dataframe(dataframe):

    index_list = list(dataframe.index.values)

    return index_list


def get_dataframe_from_multiple_index(dataframe, indeces, predictors, dv):
    dataframe_data = {}

    for idx, index_value in enumerate(indeces):
        current_index_data = dataframe.loc[dataframe.index ==
                                           index_value].values[0]
        dataframe_data[idx] = current_index_data

    current_iteration_data = pd.DataFrame.from_dict(
        data=dataframe_data, orient='index', columns=predictors + [dv])

    return current_iteration_data


def display_statistics_adjusted(model_class):
    nagelkerke_adjusted_bootstrap_ci = [
        model_class.base_ci.confidence_interval[0][0], model_class.base_ci.confidence_interval[1][0]]
    nagelkerke_adjusted = model_class.adjusted_nagelkerke

    print(f"delta R2 = {round(nagelkerke_adjusted, 4)}"
          f"[{round(nagelkerke_adjusted_bootstrap_ci[0], 4)} - {round(nagelkerke_adjusted_bootstrap_ci[1], 4)}]")

    c_stat_adjusted_bootstrap_ci = [
        model_class.base_ci.confidence_interval[0][1], model_class.base_ci.confidence_interval[1][1]]
    c_stat_adjusted = model_class.adjusted_c_stat

    print(f"delta C-statistic = {round(c_stat_adjusted, 4)}"
          f"[{round(c_stat_adjusted_bootstrap_ci[0], 4)} - {round(c_stat_adjusted_bootstrap_ci[1], 4)}]")


def create_hist_boot_distribution(data, value):

    fig, ax = plt.subplots()
    ax.hist(data, bins=25, density=True)
    # ax.plot(x, pdf)
    ax.set_title(
        f'Normal Approximation of the Bootstrap Distribution for {value}')
    ax.set_xlabel('statistic value')
    ax.set_ylabel('pdf')
    plt.show()


# All functions and class for logistic regression calculations

# given a list of predictors and interactions return a string compatible
# with statsmodels formula api

def get_predictor_string(predictors, interactions={}):

    if predictors[0] in interactions:
        predictor_string = f"{predictors[0]}:{interactions[predictors[0]]} + {predictors[0]} + "
    else:
        predictor_string = f"{predictors[0]} + "

    for predictor in predictors[1:-1]:
        if predictor in interactions:
            predictor_string += f"{predictor}:{interactions[predictor]} + {predictor} + "
        else:
            predictor_string += f"{predictor} + "

    predictor_string += predictors[-1]

    return predictor_string


def calculate_nagelkerke(y, y_pred):

    y_base = sum(y) / float(y.shape[0]) * np.ones(y.shape[0])

    log_loss_normal = -log_loss(y, y_pred)

    log_loss_base = -log_loss(y, y_base)

    cox_r2 = 1 - np.exp(((log_loss_normal) - (log_loss_base)) / len(y))
    nagelkerke = (cox_r2 / (1 - np.exp(-1 * (log_loss_normal / len(y)))))

    return nagelkerke


def calculate_nagelkerke_new(y, y_pred_probs):
    eps = 1e-15
    y_pred_probs = np.clip(y_pred_probs, eps, 1 - eps)

    # Log-likelihoods
    ll_model = np.sum(y * np.log(y_pred_probs) +
                      (1 - y) * np.log(1 - y_pred_probs))
    p_null = np.mean(y)
    ll_null = np.sum(y * np.log(p_null) +
                     (1 - y) * np.log(1 - p_null))

    n = len(y)

    # Correct Cox & Snell R²
    r2_cs = 1 - np.exp((2 * (ll_null - ll_model)) / n)

    # Maximum possible Cox & Snell R²
    r2_max = 1 - np.exp((2 * ll_null) / n)

    # Nagelkerke R²
    r2_nagelkerke = r2_cs / r2_max

    return r2_nagelkerke


def testMulticolinearity(dataframe, columns):
    testing_data = dataframe[columns].dropna(axis=0)
    vif_data = pd.DataFrame()
    vif_data["Feature"] = testing_data.columns
    vif_data["VIF"] = [variance_inflation_factor(
        testing_data.values, i) for i in range(len(testing_data.columns))]
    print(vif_data)


def likelihood_ratio_test(deviance_baseline, deviance_new, df_baseline, df_new, type):
    likelihood_ratio = deviance_baseline - deviance_new
    dof = df_new - df_baseline
    p_val = stats.chi2.sf(likelihood_ratio, dof)

    print(f"--> Likelihood ratio test {type} P-value: {round(p_val, 3)}")

    return p_val


def createRestrictedCubicSplines(dataframe, dv, predictor, display_transformed=True, create_spline_graph=True, verbose=False):
    spline_data = dataframe.loc[(dataframe[dv] != 99) & (
        dataframe[predictor].notna())]

    spline_data_transformed = dmatrix(
        "cr(train, df=4)", {"train": spline_data[predictor]}, return_type='dataframe')
    spline_data_transformed = pd.concat(
        [spline_data_transformed, dataframe], axis=1)

    # spline_data_transformed.rename({'Intercept': 'intercept', 'cr(train, df=5)[0]':f"spline_1_{predictor}", 'cr(train, df=5)[1]':f"spline_2_{predictor}", 'cr(train, df=5)[2]':f"spline_3_{predictor}", 'cr(train, df=5)[3]':f"spline_4_{predictor}", 'cr(train, df=5)[4]':f"spline_5_{predictor}"}, axis=1, inplace=True)

    spline_data_transformed.rename({'Intercept': 'intercept', 'cr(train, df=4)[0]': f"spline_1_{predictor}", 'cr(train, df=4)[1]':
                                    f"spline_2_{predictor}", 'cr(train, df=4)[2]': f"spline_3_{predictor}", 'cr(train, df=4)[3]': f"spline_4_{predictor}"}, axis=1, inplace=True)

    if display_transformed:
        spline_data_transformed.head()

    # base_model, new_model = compareModels([predictor], [f"spline_1_{predictor}",f"spline_2_{predictor}", f"spline_3_{predictor}", f"spline_4_{predictor}", f"spline_5_{predictor}"], dv,  spline_data_transformed, False, spline_data_transformed)
    base_model, new_model, lr_pvalue = compareModels([predictor], [f"spline_1_{predictor}", f"spline_2_{predictor}",
                                                                   f"spline_3_{predictor}", f"spline_4_{predictor}"], dv,  spline_data_transformed, False, spline_data_transformed, True, verbose)

    if create_spline_graph:
        xp = np.linspace(spline_data_transformed[predictor].min(
        ), spline_data_transformed[predictor].max(), spline_data_transformed[predictor].count())
        transformed_random_data = dmatrix(
            'cr(xp, df=3)', {"xp": xp}, return_type='dataframe')
        transformed_random_data.rename({'Intercept': 'intercept', 'cr(train, df=3)[0]': f"spline_1_{predictor}", 'cr(train, df=3)[1]':
                                        f"spline_2_{predictor}", 'cr(train, df=3)[2]': f"spline_3_{predictor}"}, axis=1, inplace=True)
        pred = new_model.model.predict(transformed_random_data)
        pred_basic = base_model.model.predict(xp)

        sns.scatterplot(x=xp, y=spline_data_transformed['gose_dichotimized'])
        plt.plot(
            xp, pred, label="restricted cubic spline with (3 knots)", color='orange')
        plt.plot(xp, pred_basic, color='green')
        plt.legend()

    return (spline_data_transformed, lr_pvalue)
