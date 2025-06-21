from Bio import Entrez
from dotenv import load_dotenv
import os
import sqlite3



def fetch_abstracts(disease, max_articles=400) -> list:
    handle = Entrez.esearch(db="pubmed", term=disease, retmax=max_articles)
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]

    abstracts = []
    for pmid in ids:
        summary_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstracts.append(summary_handle.read())
    return abstracts


class DataInjection:
    