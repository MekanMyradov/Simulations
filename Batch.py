import numpy as np
import matplotlib.pyplot as plt

class Batch:
    def __init__(self, x_mu: np.ndarray, x_sd: np.ndarray, y_mu: np.ndarray, y_sd: np.ndarray, weights: np.ndarray, colors: np.ndarray, proj: np.ndarray, batch_id: str, n_cells: int=1000, n_genes: int=100, gene_specific: bool=False) -> None:
        """
        Parameters:
            x_mu: x-coordinates of the components' mean.
            x_sd: standard deviation of the components.
            y_mu: y-coordinates of the components' mean.
            y_sd: standard deviation of the components.
            weights: weights of components.
            colors: colors of components.
            proj: random projection matrix; sample from standard normal distribution to project data to n_genes dimensional space.
            batch_id: id of the batch.
            n_cells: number of cells (i.e, observations or points)
            n_genes: number of genes (i.e, features)
            gene_specific: whether to add gene specific noise (batch effect) or not
        """

        self.x_mu = x_mu
        self.x_sd = x_sd

        self.y_mu = y_mu
        self.y_sd = y_sd

        self.weights = weights
        self.colors = colors
        self.proj = proj

        self.batch_id = batch_id

        self.n_cells = n_cells
        self.n_genes = n_genes

        self.gene_specific = gene_specific


    def generate_data(self):
        """
        Generate n_cells x 2 dimensional data
        """
        n_components = self.weights.size    # number of codmponents (i.e, cell types)
        component_id = np.arange(0, n_components)
        self.component_id = component_id

        # Generate a random sample with replacement
        components = np.random.choice(component_id, size=self.n_cells, replace=True, p=self.weights)    # shows which point is sampled from which component
        self.components = components

        """
        print("# 0s:", components[components == 0].size)
        print("# 1s:", components[components == 1].size)
        print("# 2s:", components[components == 2].size)
        """

        x_means = [self.x_mu[idx] for idx in components]
        y_means = [self.y_mu[idx] for idx in components]

        x_sd = [self.x_sd[idx] for idx in components]
        y_sd = [self.y_sd[idx] for idx in components]

        # Sampling locations for cells in each component
        samples = np.column_stack([
            np.random.normal(loc=x_means, scale=x_sd, size=self.n_cells),   # generate x-coordinates of cells using Gaussian distr
            np.random.normal(loc=y_means, scale=y_sd, size=self.n_cells)    # generate y-coordinates of cells using Gaussian distr
        ])

        self.samples = samples


    def save_clusters(self):
        """
        Save plot of n_cells x 2 data
        """
        plt.figure()  # Create a new figure for each plot to isolate legends

        # Create a scatter plot for each cluster
        for cluster_label in self.component_id:
            cluster_indices = (self.components == cluster_label)
            
            plt.scatter(self.samples[cluster_indices, 0], self.samples[cluster_indices, 1], c=self.colors[cluster_label], label=f'p={self.weights[cluster_label]}', marker='o', s=5)
    	
        plt.title(self.batch_id)

        plt.legend()
        plt.savefig("{}.png".format(self.batch_id))


    def project_data(self):
        """
        Project data from n_cells x 2 to n_cells x n_genes dimension
        """
        
        # (1000, 2) @ (100, 2).T = (1000, 100)
        A = np.dot(self.samples, self.proj.T)
        
        # Add standard normal noise
        noise = np.random.normal(size=(self.n_cells, self.n_genes))

        A += noise

        if self.gene_specific:
            # Add some noise (batch effect) to each gene. Noise of particular gene across cells should be the same
            gene_noise = np.random.normal(size=self.n_genes)
            
            A += gene_noise     # broadcasting

        return A
