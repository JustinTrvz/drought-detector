# Drought detection and prediction
This project is about detecting and predicting drought. The drought index used in this project is the Standardized Precipitation Evapotranspiration Index (SPEI).

# Drought detection
For the detection of a drought the SPEI is used. The SPEI indicates a drought, normal or wet condition when a specific value is reached.

| SPEI value | Condition     |
|------------|---------------|
| x <= -2.0 | Extreme drought      |
| -1.99 <= x <= -1.50  | Severe drought |
| -1.49 <= x <= -1.0  | Moderate drought |
| -0.99 <= x < 0.99   | Normal       |
| 1.0 <= x <= 1.49    | Moderately Wet          |
| 1.50 <= x <= 1.99    | Severely wet      |
| 2.0 <= x    | Extremely wet |


# Drought prediction
For the prediction of the SPEI, a Gradient Boosting Regressor is used.