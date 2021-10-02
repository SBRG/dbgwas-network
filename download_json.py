#!/usr/bin/env python3

from selenium import webdriver
from selenium.common.exceptions import ElementClickInterceptedException
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

import time
import os
import glob
import tqdm


def getjsonfile(dbgwas_file, download_file, dest_file,  chrome_loc, driver):
    """
    Opens the interactive html output of DBGWAS and downlaods the network json file to
    the file directory. Requires chrome driver to work and only works with chrome.
    """
    # assumes default download destination in downloads folder
    if os.path.isfile(download_file):
        msg = f'{download_file} file already exists. Either move or delete the file'
        raise FileExistsError(msg)
    try:
        driver.get(dbgwas_file)
        download_button = driver.find_element_by_class_name('btn-default').click()

        download_link = WebDriverWait(driver, 20).until(
                  EC.element_to_be_clickable((By.ID, "graph_cytoscape")))
        download_link.click()

        start = time.time()
        while not os.path.isfile(download_file):
            if time.time() - start > 10:
                break
            time.sleep(0.5)

        os.rename(download_file, dest_file)
    except ElementClickInterceptedException as e:
        print(e)
        print(f'Download failed for {dbgwas_file}. Download it manually instead.')
        
if __name__ == '__main__':
    import argparse
    
    ddir_help = 'path to the directory containing the component specific visualization file.\ Directory should end with visualisations/components/'
    dsdir_help = 'path to directory where all the output json files will be saved'   
        
    p = argparse.ArgumentParser(description='Download json network files from interactive DB\
                                             GWAS HTML output.')
    p.add_argument('dbgwas_dir', help=ddrir_help)
    p.add_argument('dest_dir', help=dsdir_help)
    p.add_argument('chome_loc', help='path to chrome driver')
    p.add_argument('download_loc', help='download destination of chrome profile. Default chrome\
    setting is /home/usr/downloads/')
    params = vars(p.parse_args())
    
    dbgwas_dir = os.path.abspath(params['dbgwas_dir'])
    dest_dir = os.path.abspath(params['dest_dir'])
    chrome_loc = os.path.abspath(params['chrome_loc'])
    download_loc = os.path.abspath(params['download_loc'])
    
    download_file = os.path.join(download_loc, 'graph2cytoscapeDesktop.json')
    if os.path.isfile(download_file):
        msg = f'{download_file} file already exists. Either move or delete the file'
        raise FileExistsError(msg)
    driver = webdriver.Chrome(chrome_loc)
    
    # download all the network json files
    for html_file in tqdm.tqdm(glob.glob(dbgwas_dir + '/*.html')):
        dbgwas_file = 'file://' + html_file
        comp = os.path.basename(html_file).replace('.html', '.json')
        dest_file = os.path.join(dest_dir, comp)
        if os.path.isfile(dest_file):
            continue
    #     print(f'createing {comp}')
        getjsonfile(dbgwas_file,download_file, dest_file, chrome_loc, driver)
    driver.quit()