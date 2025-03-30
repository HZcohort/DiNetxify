from DiseaseNetPy.visualization import ThreeDimensionalNetwork
import pandas as pd

if __name__ =="__main__":
    phewas_result = pd.read_csv(
        "D:/GitHubWarehouse/DiseaseNet/src/result/phewas_result.csv"
    )
    comorbidity_result = pd.read_csv(
        "D:/GitHubWarehouse/DiseaseNet/src/result/comorbidity_result.csv"
    )
    trajectory_result = pd.read_csv(
        "D:/GitHubWarehouse/DiseaseNet/src/result/trajectory_result.csv"
    )

    network = ThreeDimensionalNetwork(
        phewas_result=phewas_result,
        comorbidity_result=comorbidity_result,
        trajectory_result=trajectory_result
    )

    network.threeDimension_plot(
        "D:/GitHubWarehouse/DiseaseNet/src/figure/threeDimensionPlot"
    )
    network.significant_trajectory_plot(
        "D:/GitHubWarehouse/DiseaseNet/src/figure"
    )
    network.phewas_plot(
        "D:/GitHubWarehouse/DiseaseNet/src/figure/phewas_plot"
    )
    network.comorbidity_network_plot(
        "D:/GitHubWarehouse/DiseaseNet/src/figure/comorbidity_network_plot"
    )