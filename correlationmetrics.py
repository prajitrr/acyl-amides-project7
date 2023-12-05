#Defines correlation metrics to return correlation between two numpy arrays
#Metrics specified as lambda functions to allow for seaborn compatibility
#Contains Pearson, square root Pearson, and Spearman metrics
#Contact Prajit Rajkumar (prajkumar@ucsd.edu) for any questions

#Note that you may have to install some of the imports using 
#pip or another installer
import numpy as np

#Pearson Correlation Metric
pearson = lambda u, v: (
    np.abs(
        1 - (
            np.dot(
                np.subtract(u, np.average(u)),
                np.subtract(v, np.average(v))
            )
            /
            np.sqrt(
                (
                    (np.square(np.subtract(u, np.average(u)))).sum()
                )
                *
                (
                    (np.square(np.subtract(v, np.average(v)))).sum()
                )
            )
        )
    )
)

#Square Root Pearson Correlation Metric
sqrtpearson = lambda u, v: (
    np.sqrt(
        np.abs(
            1 - (
                np.dot(
                    np.subtract(u, np.average(u)),
                    np.subtract(v, np.average(v))
                )
                /
                np.sqrt(
                    (
                        (np.square(np.subtract(u, np.average(u)))).sum()
                    )
                    *
                    (
                        (np.square(np.subtract(v, np.average(v)))).sum()
                    )
                )
            )
        )
    )
)

#Spearman Correlation Metric
spearman = lambda u, v: (
    np.abs(
        1 - np.cov(
            ((u.argsort()).argsort()),
            ((v.argsort()).argsort())
        )[0][1]
    ) /
    (
        np.std((u.argsort()).argsort()) * np.std((v.argsort()).argsort())
    )
)