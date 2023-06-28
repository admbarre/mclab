import json
import os
import mclab
import sys
import gspread
import pandas as pd

def init_df():
    gc = gspread.oauth()
    wks = gc.open("Preps")
    sheet = wks.sheet1

    values = sheet.get_all_records()
    df = pd.DataFrame(values)
    return df

# TODO: make it so we can use other filters to subset an order
def get_rxns_info(order_date):
    df = init_df()
    todays = df[df["Prep Date"] == order_date]

    rxns = {}
    for index,prep in todays.iterrows():
        name = prep["Prep"]

        # We want to skip the preps that have too low conc
        if prep["Delta conc"] < 0:
            continue
        # TODO: hard coding p85 for today
        rxns[name] = {
            "conc" : prep["Conc [ng/ul]"],
            "seq_primers" : ["p85"]
        }
    return rxns

def load_contact_info(contact_path):
    with open(contact_path) as file:
        data = json.load(file)
    return data


if __name__ == "__main__":
    if len(sys.argv) == 2:
        _,order_date = sys.argv
        rxns_info = get_rxns_info(order_date)

        # TODO: make it so we can create new users and then select from them
        contact_path = "res/contact.json"
        contact_info = load_contact_info(contact_path)

        # TODO: printing PO for ease of copying into online order
        # This should be replaced by automatic ordering through Selenium
        print(f"PO: {contact_info['po']}")
        order = mclab.McOrder(rxns_info,contact_info)
        save_dir = os.path.abspath("orders")
        filename = order.save(save_dir)
    else:
        print("Usage: python mcorder.py <order_date>")

        
