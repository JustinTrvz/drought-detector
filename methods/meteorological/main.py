from spei_meteorological import *

start = datetime(2022, 1, 1)
end = datetime(2022, 12, 1)
out_dir = "output"
equations = ["pt", "hg", "tw", "fao56pm"]
distributions = ["p3", "ln", "ll", "gev"]

for eq in equations:
    for dist in distributions:
        print(f"Equation: '{eq}', Distribution: '{dist}'")
        calc_spei_n_months(start_date=start, end_date=end,
                          pet_equation=eq, distribution=dist,
                          save_dir=out_dir,
                          calc_1m_spei=True
                          )
