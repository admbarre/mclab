# mcdown
# Downloads all data from McLab

# Command line args
import sys

# Waiting for downloads
import time

# Loading credentials
import json

# Moving and unzipping files
import os
import shutil
import zipfile

# Selenium
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# GLOBAL DIRS
# TODO: these need to be generalized
seq_dir = "/Users/adrian/Documents/research/sequencing"
download_dir = "/Users/adrian/Downloads"
mclab_login = "https://www.mclab.com/login.php"

def load_credentials(credentials_path):
    with open(credentials_path) as file:
        data = json.load(file)
    return data

def login(credentials):
    username = credentials["username"]
    password = credentials["password"]

    driver = webdriver.Chrome()
    driver.get(mclab_login)

    # Identify login fields
    username_field = driver.find_element(By.ID,"username")
    password_field = driver.find_element(By.ID,"password")
    submit_button = driver.find_element(By.CLASS_NAME, "button-left")
    
    # Login
    username_field.send_keys(username)
    password_field.send_keys(password)
    submit_button.click()
    return driver

def download_all(driver):
    data_table_sel = "#center-main > table > tbody > tr > td > table > tbody > tr > td > table > tbody"
    data_table = driver.find_element(By.CSS_SELECTOR, data_table_sel)
    data_rows = data_table.find_elements(By.CSS_SELECTOR, "tr")
    headers,*rows = data_rows
    for row in rows:
        # The link is in the first table data entry
        link = row.find_element(By.TAG_NAME, "td")
        if is_downloaded(link.text):
            continue
        else:
            download_link(link)
    driver.quit()

def download_link(link):
    file_name = link.text
    downloaded_file_path = os.path.join(download_dir,file_name)

    link.click()
    check_interval = 15 # seconds
    max_timeout = 10 # intervals
    wait_time = check_interval * max_timeout
    interval = 0

    print("---"*15)
    print(file_name)
    print("---"*15)
    print(f"Download has started, will wait for {wait_time} seconds")
    while not os.path.exists(downloaded_file_path) and interval < max_timeout:
        time.sleep(check_interval)
        interval = interval + 1
        print(f"Waiting for download...[{interval*check_interval}/{wait_time}s]")
    # TODO: need to add check to make sure that file downloaded
    print("***")
    move_download(downloaded_file_path)

def is_downloaded(datafile):
    folder_name = datafile.split(".zip")[0]
    seqs = [f for f in os.listdir(seq_dir)]
    if folder_name in seqs:
        print(f"{datafile} already in sequencing dir")
        return True
    return False

def goto_data_page(driver):
    max_pageload = 10 #seconds
    # The selector for the download page
    download_sel = "#header > div > div.tabs > ul > li.first > a"

    # Wait for download sequences link to load
    wait = WebDriverWait(driver,max_pageload)
    # TODO: what does element do here? does it return the element?
    element = wait.until(EC.visibility_of_element_located((By.CSS_SELECTOR,download_sel)))
    download_link = driver.find_element(By.CSS_SELECTOR,download_sel)
    download_link.click()

# Does not check for file existing because precondition
def move_download(downloaded_file_path):
    filename = downloaded_file_path.split("/")[-1]
    new_path_zip = f"{seq_dir}/{filename}"
    new_dir = new_path_zip.split(".zip")[0]

    shutil.move(downloaded_file_path,new_path_zip)
    with zipfile.ZipFile(new_path_zip, "r") as zip_ref:
        zip_ref.extractall(new_dir)
    os.remove(new_path_zip)

# TODO: pick a better name for this
def main(credentials_path):
    creds = load_credentials(credentials_path)
    driver = login(creds)
    goto_data_page(driver)
    download_all(driver)

if __name__ == "__main__":
    if len(sys.argv) ==2:
        _,credentials_path = sys.argv
        main(credentials_path)
    else:
        print("Usage: python mcdown.py <mclabs credentials json>")
