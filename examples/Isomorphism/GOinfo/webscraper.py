#!/usr/bin/env python3
# coding: utf-8
"""Scrape every Physics PhD since 1970."""

from selenium import webdriver
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
# from selenium.webdriver.support.ui import WebDriverWait
# from selenium.webdriver.support import expected_conditions as EC
# from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException, TimeoutException
import time

caps = DesiredCapabilities.FIREFOX
caps["marionette"] = True

caps["binary"] = "/usr/bin/firefox"
browser = webdriver.Firefox(capabilities=caps)
time.sleep(0.5)

with open("GOterms.dat") as flines:
    GOterms = [line.strip().split(" ") for line in flines][1:]


GOtermWithAnnotation = []
try:
    for GOterm in GOterms:

        url = "http://amigo.geneontology.org/amigo/term/{}".format(GOterm[1])
        browser.get(url)
        time.sleep(0.5)
        try:

            data = browser.find_element_by_xpath("/html/body/div[2]/div[2]/div[2]/dl/dd[2]").text
        except NoSuchElementException:
            continue
        except TimeoutException:
            continue
        GOterm.append(data)
        GOtermWithAnnotation.append(GOterm)

    with open("GOtermswithannotation.dat", mode='w') as flines:
        flines.write("\n".join(" ".join(x)) for x in GOtermWithAnnotation)

except KeyboardInterrupt:
    with open("GOtermswithannotation.dat", mode='w') as flines:
        flines.write("\n".join([" ".join(x) for x in GOtermWithAnnotation]))
