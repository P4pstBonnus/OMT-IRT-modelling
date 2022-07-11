import datetime as dt
import pandas as pd

base_df_relevant = pd.read_csv("original_omt_data_motive.csv", sep=",")

final_df = pd.DataFrame(
    {"Person": [], "Picture": [], "Code": [], 
     "DAFF": [], "DACH": [], "DPOW": [], #"D0": [], 
     "TDAFF": [], "TDACH": [], "TDPOW": [], #"TD0": [],
     "SDAFF": [], "SDACH": [], "SDPOW": [], #"SD0": [],
     "Response1": [], "Response2": []
    }
)
possible_motives = ["AFF", "ACH", "POW", "0"]
picture_columns = [col for col in list(base_df_relevant.columns) 
                   if col not in ["id", "text_len_total"]]
for rowcount, row in enumerate(base_df_relevant.iterrows()):
    if rowcount % 1000 == 0:
        print(f"{dt.datetime.now()}: {rowcount}")
    df_per_person = pd.DataFrame()
    last_motive = None
    aff_times = 0
    ach_times = 0
    pow_times = 0
    zero_times = 0
    for col in base_df_relevant[picture_columns].columns:
        current_motive = row[1][col]
        results = get_vectors_and_responses_for_options_and_winner(possible_motives, current_motive)
        for i in range(3):
            df_per_person = df_per_person.append(
                pd.DataFrame({
                    "Person": [int(row[1]["id"])],
                    "Picture": [col],
                    "Code": [current_motive],
                    "DAFF": [results[i][0][0]],
                    "DACH": [results[i][0][1]],
                    "DPOW": [results[i][0][2]],
                    "D0": [results[i][0][3]],
                    "TDAFF": [1 if last_motive == "AFF" else 0],
                    "TDACH": [1 if last_motive == "ACH" else 0],
                    "TDPOW": [1 if last_motive == "POW" else 0],
                    #"TD0": [],
                    "SDAFF": [aff_times],
                    "SDACH": [ach_times],
                    "SDPOW": [pow_times],
                    "SD0": [zero_times],
                    "Response1": [results[i][1][0]],
                    "Response2": [results[i][1][1]]
                }),
                ignore_index=True
            )
        aff_times += 1 if current_motive == "AFF" else 0
        ach_times += 1 if current_motive == "ACH" else 0
        pow_times += 1 if current_motive == "POW" else 0
        zero_times += 1 if current_motive == "0" else 0
        last_motive = current_motive
    final_df = final_df.append(df_per_person, ignore_index=True)
final_df = final_df.astype(
    {"Person": int, "DAFF": int, "DACH": int, "DPOW": int, "D0": int, 
     "TDAFF": int, "TDACH": int, "TDPOW": int,
     "SDAFF": int, "SDACH": int, "SDPOW": int, "SD0": int,
     "Response1": int, "Response2": int
    }
)
