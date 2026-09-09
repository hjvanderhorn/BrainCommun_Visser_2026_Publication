from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import KernelPCA
from factor_analyzer.rotator import Rotator
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Tuple, List


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by mixing it with white.
    Input can be matplotlib color string, hex string, or RGB tuple.
    amount=0 returns original color, amount=1 returns white.
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = np.array(mc.to_rgb(c))
    white = np.array([1, 1, 1])
    return tuple((1 - amount) * c + amount * white)


class create_PCA:
    def __init__(self, features, outcome_variable, dataframe, sample_id):
        self.original_data = dataframe
        self.predictors = features
        self.outcome_variable = outcome_variable
        self.data = dataframe[features +
                              [outcome_variable, sample_id]].dropna(axis=0)
        self.features = self.data[features]
        self.labels = self.data[outcome_variable]
        self.sample_label = self.data[sample_id]

    def normalize_features(self):
        self.scaler = StandardScaler()
        self.features_scaled = self.scaler.fit_transform(self.features)

    def fit_pca(self, n_components, linear=True, kernel='rbf', gamma=0.1):
        self.n_components = n_components
        if not hasattr(self, 'features_scaled') or self.features_scaled is None:
            self.normalize_features()
        if linear:
            self.pca = PCA(n_components=self.n_components)

        if not linear:
            self.pca = KernelPCA(
                n_components=self.n_components, kernel=kernel, gamma=gamma)

        self.pc = self.pca.fit_transform(self.features_scaled)
        component_names = [f'PC{i+1}' for i in range(self.n_components)]
        self.pca_df = pd.DataFrame(data=self.pc, columns=component_names)
        self.pca_df['ct_injury_classification'] = self.labels.reset_index(
            drop=True)
        self.pca_df['record_id'] = self.sample_label.reset_index(drop=True)
        # display(self.pca_df)
        return (self.pca_df)

    def merge_original_data(self):
        merged_data = pd.merge(left=self.original_data, right=self.pca_df,
                               left_on='record_id', right_on='record_id', how='left')
        return merged_data

    def create_heatmap(self, n_components: int = 3, figsize: Tuple[int, int] = (6, 10), dpi: int = 300):
        """
        Create a heatmap of PCA loading scores.

        Parameters:
        -----------
        n_components : int, optional
            Number of principal components to include in the heatmap (default is 3)
        figsize : tuple of int, optional
            Figure size in inches (width, height) (default is (6, 10))
        dpi : int, optional
            Figure resolution (default is 300)
        """
        loading_scores = self.pca.components_[:n_components]
        loading_scores_df = pd.DataFrame(
            loading_scores,
            columns=self.features.columns,
            index=[f'PC{i+1}' for i in range(n_components)]
        )

        # Create custom viridis-based colormap for diverging data
        # Negative values (yellow) to positive values (green)
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.pyplot as plt

        # Define colors: yellow for negative, green for positive
        colors = ['#440154', '#31688e', '#35b779', '#fde725']  # viridis colors
        n_bins = 256
        cmap = LinearSegmentedColormap.from_list(
            'viridis_div', colors, N=n_bins)

        plt.figure(figsize=figsize, dpi=dpi)
        sns.heatmap(loading_scores_df.T, annot=True, cmap=cmap,
                    center=0, cbar_kws={'label': 'Loading Score'})
        plt.title('PCA Loading Scores Heatmap')
        plt.xlabel('Principal Components')
        plt.ylabel('Features')
        plt.show()

    def bootstrap_loadings(self, n_iterations):
        loading_distributions = []

        for _ in range(n_iterations):
            resampled_data = self.features.sample(
                n=len(self.features), replace=True)
            scaler = StandardScaler()
            resampled_data_scaled = scaler.fit_transform(resampled_data)
            pca = PCA(n_components=self.n_components)
            pca.fit(resampled_data_scaled)
            loading_distributions.append(pca.components_)

        loading_distributions = np.array(loading_distributions)
        self.average_loadings = np.mean(loading_distributions, axis=0)

    def bootstrap_loadings_varimax(self, n_iterations):
        bootrap_loading_distributions = []

        for _ in range(n_iterations):
            resampled_data = self.features.sample(
                n=len(self.features), replace=True)
            scaler = StandardScaler()
            resampled_data_scaled = scaler.fit_transform(resampled_data)
            pca = PCA(n_components=self.n_components)
            pca.fit(resampled_data_scaled)
            loadings = pca.components_.T
            rotator = Rotator(method='varimax')
            rotated_coeff = rotator.fit_transform(loadings)
            bootrap_loading_distributions.append(rotated_coeff)

        bootrap_loading_distributions = np.array(bootrap_loading_distributions)
        self.average_loadings_varimax = np.mean(
            bootrap_loading_distributions, axis=0)

    def create_average_loadings_heatmap(self, bootstrap=False):

        if bootstrap:
            loadings = self.average_loadings_varimax.T
        else:
            loadings = self.average_loadings
        self.average_loadings_df = pd.DataFrame(loadings, columns=self.features.columns, index=[
                                                f'PC{i+1}' for i in range(len(loadings))])
        cmap = sns.diverging_palette(220, 20, as_cmap=True)
        plt.figure(figsize=(6, 10), dpi=300)
        sns.heatmap(self.average_loadings_df.T, annot=True, cmap=cmap,
                    center=0, cbar_kws={'label': 'Average Loading Score'})
        plt.title('Average PCA Loading Scores Heatmap')
        plt.xlabel('Principal Components')
        plt.ylabel('Features')
        plt.show()

    def create_scatter_plot(self):
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=self.pca_df, x='PC1', y='PC2',
                        hue='ct_injury_classification', palette=['blue', 'red'])
        # plt.ylim((0, 10))
        # plt.xlim((-4, -3))
        plt.title('PCA of Biomarker Profiles')
        plt.xlabel(f'PC1')
        plt.ylabel(f'PC2 ')
        plt.legend(title='ct_injury_classification')
        plt.grid(True)
        plt.show()

    def varimax_rotation(self, features):
        coeff = self.pca.components_.T
        score = self.pca.transform(self.features_scaled)
        rotator = Rotator(method='varimax')
        rotated_coeff = rotator.fit_transform(coeff[:, :3])
        self.rotation_matrix = rotator.rotation_
        rotated_score = np.dot(score[:, :3], self.rotation_matrix)
        rotated_coeff_df = pd.DataFrame(rotated_coeff, index=features, columns=[
                                        'rotated_pc1', 'rotated_pc2', 'rotated_pc3'])
        rotated_score_df = pd.DataFrame(
            rotated_score, columns=['rotated_pc1', 'rotated_pc2', 'rotated_pc3'])
        return rotated_coeff_df, rotated_score_df

    def create_scree_plot(self):
        explained_variance_ratio = self.pca.explained_variance_ratio_
        plt.figure(figsize=(10, 6))
        plt.plot(range(1, len(explained_variance_ratio) + 1),
                 explained_variance_ratio, 'o-', linewidth=2, color='blue')
        plt.title('Scree Plot')
        plt.xlabel('Principal Component')
        plt.ylabel('Explained Variance Ratio')
        plt.xticks(range(1, len(explained_variance_ratio) + 1))
        plt.grid()
        plt.show()

    def create_kernel_scree_plot(self):
        # For Kernel PCA, we use eigenvalues_
        eigenvalues = self.pca.eigenvalues_
        # Normalize the eigenvalues
        explained_variance_ratio = eigenvalues / eigenvalues.sum()

        plt.figure(figsize=(10, 6))
        plt.plot(range(1, len(explained_variance_ratio) + 1),
                 explained_variance_ratio, 'o-', linewidth=2, color='blue')
        plt.title('Kernel PCA Scree Plot')
        plt.xlabel('Principal Component')
        plt.ylabel('Normalized Eigenvalue')
        plt.xticks(range(1, len(explained_variance_ratio) + 1))
        plt.grid()
        plt.show()

    def group_difference(self, components=["PC1"], rotated=False):

        if (rotated):
            self.outcome_df = self.rotated_score_df.loc[
                self.rotated_score_df[self.outcome_variable] != 99]
        else:
            self.outcome_df = self.pca_df.loc[self.pca_df[self.outcome_variable] != 99]

        print(self.outcome_variable)

        for component in components:
            shapiro_result = stats.shapiro(self.outcome_df[component])
            mannwhitneyu = stats.mannwhitneyu(x=self.outcome_df.loc[self.outcome_df[self.outcome_variable] == "COVID19", [
                component]], y=self.outcome_df.loc[self.outcome_df[self.outcome_variable] == "Control", [component]])
            print(
                f'Assessing group differences and distribution for {component}\nShapiro: {shapiro_result.pvalue}\n(mannWhitney): {mannwhitneyu.pvalue}\n')

    def group_difference_barchart(self, component):
        plt.figure(figsize=(8, 6))
        sns.barplot(x='Group', y=component,
                    data=self.outcome_df, capsize=.2)
        plt.title('PC1 loading scores by group')
        plt.xlabel('GOSE')
        plt.ylabel(f'{component} score')
        plt.axhline(0, color='gray', linewidth=0.8)
        # plt.ylim(-1, 1)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.show()

    def perform_pca_single_timepoint(
        self,
        biomarkers: List[str],
        n_components: int = 3
    ) -> pd.DataFrame:
        """
        Performs PCA analysis on a single timepoint dataframe.

        Parameters:
        -----------
        biomarkers: List[str]
            List of biomarker names
        n_components: int
            Number of components to use

        Returns:
        --------
        pd.DataFrame
            DataFrame containing the PCA loadings
        """

        total_components = n_components
        self.fit_pca(total_components)

        component_names = [f'PC{i+1}' for i in range(total_components)]
        loadings = self.pca.components_.T
        loadings_df = pd.DataFrame(
            loadings, index=biomarkers, columns=component_names)

        rotated_coeff_df, self.rotated_score_df = self.varimax_rotation(
            biomarkers)

        label_columns = ['Group', 'sample_id']
        acute_labels = self.data[label_columns]
        self.rotated_score_df[label_columns] = acute_labels.reset_index(
            drop=True)
        self.rotated_score_df['rotated_pc2'] *= -1
        self.rotated_score_df['rotated_pc3'] *= -1

        return loadings_df, self.rotated_score_df, rotated_coeff_df

    def perform_pca_analysis(
        self,
        acute_data: pd.DataFrame,
        subacute_data: pd.DataFrame,
        biomarkers: List[str],
        hc_acute_data: pd.DataFrame,
        n_components: int = 3
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Performs PCA analysis on acute and subacute data, including varimax rotation.

        Parameters:
        -----------
        acute_data : pd.DataFrame
            DataFrame containing acute timepoint data
        subacute_data : pd.DataFrame
            DataFrame containing subacute timepoint data
        biomarkers : list
            List of biomarker column names
        hc_acute_data : pd.DataFrame
            Healthy control acute data for normalization
        n_components : int, optional
            Number of components to use (default is 3)

        Returns:
        --------
        tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
            - Rotated coefficients DataFrame
            - Combined rotated scores DataFrame
            - PCA loadings DataFrame
        """
        # Store the original n_components for loadings
        total_components = n_components

        # Fit PCA with all requested components
        self.fit_pca(total_components)

        # Get loadings for all components
        component_names = [f'PC{i+1}' for i in range(total_components)]
        loadings = self.pca.components_.T
        loadings_df = pd.DataFrame(
            loadings, index=biomarkers, columns=component_names)

        # Perform varimax rotation (always on first 3 components)
        rotated_coeff_df, rotated_score_df = self.varimax_rotation(biomarkers)

        # Adjust signs for consistency
        rotated_coeff_df['rotated_pc2'] *= -1
        rotated_coeff_df['rotated_pc3'] *= -1

        # Add labels to rotated scores
        label_columns = ['gose_dichotimized',
                         'timePoint', 'record_id', 'group']
        acute_labels = acute_data[label_columns]
        rotated_score_df[label_columns] = acute_labels.reset_index(drop=True)
        rotated_score_df['rotated_pc2'] *= -1
        rotated_score_df['rotated_pc3'] *= -1

        # Process subacute data
        subacute_labels = subacute_data[label_columns]
        subacute_pca_data = subacute_data[biomarkers]

        # Normalize subacute data
        hc_acute_pca_data = hc_acute_data[biomarkers]
        # means = hc_acute_pca_data.mean(axis=0)
        # stds = hc_acute_pca_data.std(axis=0)
        # subacute_normalized = (subacute_pca_data - means) / stds
        subacute_normalized = self.scaler.transform(subacute_pca_data)

        # Fit existing PCA to subacute data (use only first 3 components for rotation)
        fitted_subacute = self.pca.transform(subacute_normalized)
        rotated_score_subacute = np.dot(
            fitted_subacute[:, :3], self.rotation_matrix)

        # Create subacute rotated scores DataFrame
        rotated_score_df_subacute = pd.DataFrame(
            rotated_score_subacute,
            columns=['rotated_pc1', 'rotated_pc2', 'rotated_pc3']
        )
        rotated_score_df_subacute['rotated_pc2'] *= -1
        rotated_score_df_subacute['rotated_pc3'] *= -1
        rotated_score_df_subacute[label_columns] = subacute_labels.reset_index(
            drop=True)

        # Combine acute and subacute data
        self.combined_rotated_df = pd.concat(
            [rotated_score_df, rotated_score_df_subacute], axis=0)
        self.combined_rotated_df.reset_index(
            drop=True, inplace=True)

        return rotated_coeff_df, self.combined_rotated_df, loadings_df

    def perform_kernel_pca(self, n_components: int = 3, kernel: str = 'rbf', gamma: float = 0.1, biomarkers: List[str] = None, acute_data: pd.DataFrame = None, subacute_data: pd.DataFrame = None):
        """
        Perform Kernel PCA, calculate feature correlations, and project both acute and subacute data.

        Parameters:
        -----------
        n_components : int
            Number of components to extract
        kernel : str
            Kernel type ('rbf', 'linear', 'poly', etc.)
        gamma : float
            Kernel coefficient for 'rbf', 'poly' and 'sigmoid' kernels
        biomarkers : List[str]
            List of biomarker names
        acute_data : pd.DataFrame
            Acute timepoint data
        subacute_data : pd.DataFrame
            Subacute timepoint data

        Returns:
        --------
        Tuple[pd.DataFrame, pd.DataFrame]
            - DataFrame containing the correlation coefficients between features and components
            - DataFrame containing the projected data (both acute and subacute) with labels
        """
        total_components = n_components

        # Fit PCA with all requested components
        self.fit_pca(total_components, linear=False,
                     kernel=kernel, gamma=gamma)

        component_names = [f'PC{i+1}' for i in range(total_components)]

        # Calculate correlations between original features and KPCA components
        correlations = []
        for i in range(total_components):
            component_corr = [stats.pearsonr(self.features[col], self.pc[:, i])[0]
                              for col in self.features.columns]
            correlations.append(component_corr)

        # Create correlation matrix DataFrame
        loadings_df = pd.DataFrame(
            correlations,
            columns=self.features.columns,
            index=component_names
        )

        # Project acute data
        acute_projected = pd.DataFrame(
            self.pc,
            columns=component_names
        )

        # Add labels for acute data
        label_columns = ['gose_dichotimized',
                         'timePoint', 'record_id', 'group']
        if acute_data is not None:
            acute_labels = acute_data[label_columns].reset_index(drop=True)
            acute_projected = pd.concat(
                [acute_projected, acute_labels], axis=1)

        # Project subacute data if provided
        if subacute_data is not None and biomarkers is not None:
            # Prepare subacute data
            subacute_features = subacute_data[biomarkers]

            # Normalize subacute data using the same scaler
            subacute_scaled = self.scaler.transform(subacute_features)

            # Project subacute data using the fitted kernel PCA
            subacute_pc = self.pca.transform(subacute_scaled)

            # Create DataFrame for subacute projected data
            subacute_projected = pd.DataFrame(
                subacute_pc,
                columns=component_names
            )

            # Add labels for subacute data
            subacute_labels = subacute_data[label_columns].reset_index(
                drop=True)
            subacute_projected = pd.concat(
                [subacute_projected, subacute_labels], axis=1)

            # Combine acute and subacute projected data
            projected_df = pd.concat(
                [acute_projected, subacute_projected], axis=0)
            projected_df.reset_index(drop=True, inplace=True)
        else:
            projected_df = acute_projected

        return loadings_df, projected_df

    def create_rotated_heatmap(self, figsize: Tuple[int, int] = (6, 10), dpi: int = 300, sort_by: str = 'Rotated PC3', label_dict: dict = None):
        """
        Create a heatmap of varimax rotated loading scores.

        Parameters:
        -----------
        figsize : tuple of int, optional
            Figure size in inches (width, height) (default is (6, 10))
        dpi : int, optional
            Figure resolution (default is 300)
        sort_by : str, optional
            Column name to use for sorting the y-axis (default is 'Rotated PC3')
        label_dict : dict, optional
            Dictionary mapping original feature names to custom labels for y-axis
        """
        if not hasattr(self, 'rotation_matrix'):
            raise AttributeError(
                "No rotation matrix found. Please run varimax_rotation() first.")

        coeff = self.pca.components_.T
        rotator = Rotator(method='varimax')
        rotated_coeff = rotator.fit_transform(coeff[:, :3])

        # Create DataFrame for the rotated coefficients
        rotated_scores_df = pd.DataFrame(
            rotated_coeff,
            index=self.features.columns,
            columns=['Rotated PC1', 'Rotated PC2', 'Rotated PC3']
        )

        # Adjust signs for consistency with the rest of the analysis
        rotated_scores_df['Rotated PC2'] *= -1
        rotated_scores_df['Rotated PC3'] *= -1

        # Sort the DataFrame by the specified column
        rotated_scores_df_sorted = rotated_scores_df.sort_values(
            by=sort_by, ascending=False)

        # Apply custom labels if provided
        if label_dict is not None:
            rotated_scores_df_sorted.index = rotated_scores_df_sorted.index.map(
                lambda x: label_dict.get(x, x))

        # Create custom diverging colormap: yellow (negative) -> white (zero) -> green (positive)
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.pyplot as plt

        # Define colors: yellow to white to green (lightened)
        colors = [lighten_color('#FFD700', 0.3), '#FFFFFF', lighten_color(
            '#228B22', 0.3)]  # yellow and green lightened
        n_bins = 256
        cmap = LinearSegmentedColormap.from_list(
            'yellow_white_green', colors, N=n_bins)

        plt.figure(figsize=figsize, dpi=dpi)
        sns.heatmap(rotated_scores_df_sorted, annot=True, cmap=cmap,
                    center=0, cbar_kws={'label': 'Rotated Loading Score'})
        plt.xlabel('Rotated Principal Components')
        plt.ylabel('RSI metabolites')
        plt.yticks(fontweight='bold')
        plt.tick_params(axis='y', which='minor', length=0)
        plt.tight_layout()
        plt.savefig(r'/users/koen_phd/documents/phd/biomarkers/rsi/rotated_heatmap.svg',
                    dpi=300, bbox_inches='tight')
        plt.show()

        return rotated_scores_df_sorted

    def create_kernel_loadings_heatmap(self, n_components: int = 3, threshold: float = 0.4, figsize: Tuple[int, int] = (6, 10), dpi: int = 300):
        """
        Create a heatmap of kernel PCA feature contributions using correlation analysis.

        Parameters:
        -----------
        n_components : int, optional
            Number of components to include in the heatmap (default is 3)
        threshold : float, optional
            Absolute correlation threshold for considering significant loadings (default is 0.4)
        figsize : tuple of int, optional
            Figure size in inches (width, height) (default is (6, 10))
        dpi : int, optional
            Figure resolution (default is 300)
        """
        # Calculate correlations between original features and KPCA components
        correlations = []
        for i in range(n_components):
            component_corr = [stats.pearsonr(self.features[col], self.pc[:, i])[0]
                              for col in self.features.columns]
            correlations.append(component_corr)

        # Create correlation matrix DataFrame
        loading_scores_df = pd.DataFrame(
            correlations,
            columns=self.features.columns,
            index=[f'PC{i+1}' for i in range(n_components)]
        )

        # Create heatmap
        cmap = sns.diverging_palette(220, 20, as_cmap=True)
        plt.figure(figsize=figsize, dpi=dpi)
        sns.heatmap(loading_scores_df.T, annot=True, cmap=cmap,
                    center=0, cbar_kws={'label': 'Correlation Coefficient'})
        plt.title('Kernel PCA Feature Correlations Heatmap')
        plt.xlabel('Principal Components')
        plt.ylabel('Features')
        plt.tight_layout()
        plt.show()

        # Print significant features for each component
        for i in range(n_components):
            significant_features = loading_scores_df.iloc[i][
                abs(loading_scores_df.iloc[i]) >= threshold
            ]
            if len(significant_features) > 0:
                print(
                    f"\nSignificant features for PC{i+1} (|correlation| >= {threshold}):")
                for feat, corr in significant_features.items():
                    print(f"{feat}: {corr:.3f}")

        return loading_scores_df

    def kernel_group_difference(self, projected_data: pd.DataFrame, components: List[str] = ["PC1"]):
        """
        Calculate group differences for specified kernel PCA components.

        Parameters:
        -----------
        projected_data : pd.DataFrame
            DataFrame containing the projected data and labels
        components : List[str]
            List of components to analyze (e.g., ["PC1"], ["PC1", "PC2"])
        """
        for component in components:
            # Remove any missing values in outcome
            component_data = projected_data.loc[
                (projected_data['gose_dichotimized'] != 99) &
                (projected_data['group'] != 0)
            ]

            # Normality test
            shapiro_result = stats.shapiro(component_data[component])

            # Group difference test
            mannwhitney = stats.mannwhitneyu(
                x=component_data.loc[component_data['gose_dichotimized']
                                     == 1, component],
                y=component_data.loc[component_data['gose_dichotimized']
                                     == 0, component],
                alternative='two-sided'
            )

            # Calculate effect size (rank-biserial correlation)
            n1 = sum(component_data['gose_dichotimized'] == 1)
            n0 = sum(component_data['gose_dichotimized'] == 0)
            effect_size = 1 - (2 * mannwhitney.statistic / (n1 * n0))

            # Calculate descriptive statistics
            desc_stats = component_data.groupby('gose_dichotimized')[component].agg([
                'count', 'mean', 'std', 'median'
            ])

            print(f"\nAnalysis for {component}:")
            print("-" * 50)
            print("Descriptive Statistics:")
            print(desc_stats)
            print("\nNormality Test (Shapiro-Wilk):")
            print(f"p-value: {shapiro_result.pvalue:.4f}")
            print("\nGroup Difference Test (Mann-Whitney U):")
            print(f"p-value: {mannwhitney.pvalue:.4f}")
            print(
                f"Effect size (rank-biserial correlation): {effect_size:.4f}")

            # Create violin plot
            plt.figure(figsize=(8, 6))
            sns.violinplot(
                data=component_data,
                x='gose_dichotimized',
                y=component,
                inner='box'
            )
            plt.title(f'{component} Distribution by Group')
            plt.xlabel('GOSE Dichotomized')
            plt.ylabel(f'{component} Score')
            plt.show()

    def compare_hc_mtbi(self, projected_data: pd.DataFrame, components: List[str] = ["PC1"]):
        """
        Compare kernel PCA components between HC (group=0) and mTBI (group=1) at both timepoints.

        Parameters:
        -----------
        projected_data : pd.DataFrame
            DataFrame containing the projected data and labels
        components : List[str]
            List of components to analyze (e.g., ["PC1"], ["PC1", "PC2"])
        """
        for component in components:
            # Remove any missing values
            component_data = projected_data.dropna(
                subset=[component, 'group', 'timePoint'])

            # Create group labels
            component_data['group_time'] = 'HC'  # Default for group 0
            # Label mTBI acute (group 1, timePoint 1)
            component_data.loc[(component_data['group'] == 1) &
                               (component_data['timePoint'] == 1), 'group_time'] = 'mTBI_acute'
            # Label mTBI subacute (group 1, timePoint 2)
            component_data.loc[(component_data['group'] == 1) &
                               (component_data['timePoint'] == 2), 'group_time'] = 'mTBI_subacute'

            # Perform Kruskal-Wallis test (for overall comparison)
            kruskal = stats.kruskal(
                component_data[component_data['group_time']
                               == 'HC'][component],
                component_data[component_data['group_time']
                               == 'mTBI_acute'][component],
                component_data[component_data['group_time']
                               == 'mTBI_subacute'][component]
            )

            # Perform pairwise Mann-Whitney U tests
            pairs = [
                ('HC', 'mTBI_acute'),
                ('HC', 'mTBI_subacute'),
                ('mTBI_acute', 'mTBI_subacute')
            ]

            pairwise_tests = {}
            for group1, group2 in pairs:
                test = stats.mannwhitneyu(
                    x=component_data[component_data['group_time']
                                     == group1][component],
                    y=component_data[component_data['group_time']
                                     == group2][component],
                    alternative='two-sided'
                )
                pairwise_tests[(group1, group2)] = test

            # Calculate descriptive statistics
            desc_stats = component_data.groupby('group_time')[component].agg([
                'count', 'mean', 'std', 'median'
            ])

            # Print results
            print(f"\nAnalysis for {component}")
            print("-" * 50)
            print("\nDescriptive Statistics:")
            print(desc_stats)
            print(f"\nKruskal-Wallis test (overall comparison):")
            print(f"H-statistic: {kruskal.statistic:.4f}")
            print(f"p-value: {kruskal.pvalue:.4f}")
            print("\nPairwise Mann-Whitney U tests:")
            for (group1, group2), test in pairwise_tests.items():
                print(f"\n{group1} vs {group2}:")
                print(f"p-value: {test.pvalue:.4f}")

            # Create violin plot
            plt.figure(figsize=(10, 6))
            sns.violinplot(
                data=component_data,
                x='group_time',
                y=component,
                inner='box',
                order=['HC', 'mTBI_acute', 'mTBI_subacute']
            )
            plt.title(f'{component} Distribution: HC vs mTBI (Acute/Subacute)')
            plt.xlabel('Group')
            plt.ylabel(f'{component} Score')
            plt.xticks(rotation=45)

            # Add statistical annotation
            if kruskal.pvalue < 0.001:
                stat_text = 'Kruskal-Wallis: p < 0.001'
            else:
                stat_text = f'Kruskal-Wallis: p = {kruskal.pvalue:.3f}'
            plt.text(1, plt.ylim()[1], stat_text,
                     horizontalalignment='center',
                     verticalalignment='bottom')

            plt.tight_layout()
            plt.show()

            # Create boxplot with individual points
            plt.figure(figsize=(10, 6))
            sns.boxplot(
                data=component_data,
                x='group_time',
                y=component,
                width=0.5,
                order=['HC', 'mTBI_acute', 'mTBI_subacute']
            )
            sns.swarmplot(
                data=component_data,
                x='group_time',
                y=component,
                color='.25',
                size=4,
                alpha=0.5,
                order=['HC', 'mTBI_acute', 'mTBI_subacute']
            )
            plt.title(f'{component} Distribution: HC vs mTBI (Acute/Subacute)')
            plt.xlabel('Group')
            plt.ylabel(f'{component} Score')
            plt.xticks(rotation=45)

            # Add statistical annotation
            plt.text(1, plt.ylim()[1], stat_text,
                     horizontalalignment='center',
                     verticalalignment='bottom')

            plt.tight_layout()
            plt.show()
