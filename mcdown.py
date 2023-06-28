# mcdown
# donwloads the latest data from McLab's website

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

def download_latest(credentials):
    max_pageload = 10 #seconds
    # The selector for the download page
    download_sel = "#header > div > div.tabs > ul > li.first > a"
    # The selector for the top most download link
    # TODO: probably a better way to read this table...?
    # we could probably download the full table and then just update based on
    # the date
    top_download_sel = "#center-main > table > tbody > tr > td > table > tbody > tr > td > table > tbody > tr:nth-child(2) > td:nth-child(1) > a"

    driver = login(credentials)
    # Wait for download sequences link to load
    wait = WebDriverWait(driver,max_pageload)
    element = wait.until(EC.visibility_of_element_located((By.CSS_SELECTOR,download_sel)))
    download_link = driver.find_element(By.CSS_SELECTOR,download_sel)
    download_link.click()

    # Wait for latest data link to load
    wait = WebDriverWait(driver,max_pageload)
    element = wait.until(EC.visibility_of_element_located((By.CSS_SELECTOR,top_download_sel)))
    download_link = driver.find_element(By.CSS_SELECTOR,top_download_sel)
    file_name = download_link.text
    file_path = os.path.join(download_dir,file_name)
    new_path = f"{seq_dir}/{file_name}"
    new_dir = new_path.split(".zip")[0]

    # TODO: there is a lot of repetition here and room for optimizing the flow
    if os.path.exists(file_path):
        print("Already downloaded latest")
        driver.quit()
        return file_path
    elif os.path.exists(new_dir):
        print("Already downloaded and moved latest")
        driver.quit()
        return file_path

    download_link.click()
    check_interval = 15 # seconds
    max_timeout = 10 # intervals
    wait_time = check_interval * max_timeout
    interval = 0

    print(f"Download has started, will wait for {wait_time} seconds")
    while not os.path.exists(file_path) and interval < max_timeout:
        time.sleep(check_interval)
        interval = interval + 1
        print(f"Waiting for download...[{interval*check_interval}/{wait_time}s]")
    driver.quit()
    return file_path

def move_data(filepath):
    filename = filepath.split("/")[-1]
    new_path = f"{seq_dir}/{filename}"
    new_dir = new_path.split(".zip")[0]

    if os.path.exists(new_dir):
        print("Already moved file")
        return new_dir

    shutil.move(filepath,new_path)
    with zipfile.ZipFile(new_path, "r") as zip_ref:
        zip_ref.extractall(new_dir)
    return new_dir

# TODO: pick a better name for this
def main(credentials_path):
    credentials = load_credentials(credentials_path)
    filepath = download_latest(credentials)
    reads_dir = move_data(filepath)
    return reads_dir

if __name__ == "__main__":
    if len(sys.argv) ==2:
        _,credentials_path = sys.argv
        main(credentials_path)
    else:
        print("Usage: python mcdown.py <mclabs credentials json>")
