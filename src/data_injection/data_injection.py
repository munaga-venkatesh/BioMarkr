import os
import re
import sys
import sqlite3
from tqdm import tqdm
from typing import List
from Bio import Entrez
from dotenv import load_dotenv

from src.utils.logger import logging
from src.utils.exception import CustomException

load_dotenv()


class DataInjection:
    """
    A class to interact with PubMed API and manage abstracts
    related to a specific disease by saving them into a local SQLite database.
    """

    def __init__(self) -> None:
        """
        Initializes Entrez API credentials and sets up data directory paths.
        """
        self.data_folder_path = os.path.join(os.getcwd(), "data")
        self.db_path = os.path.join(self.data_folder_path, "abstracts.db")
        self.email = os.getenv("ENTREZ_EMAIL")
        self.apikey = os.getenv("ENTREZ_API_KEY")

        if not self.email or not self.apikey:
            raise ValueError("ENTREZ_EMAIL or ENTREZ_API_KEY not found in environment variables.")

        Entrez.email = self.email
        Entrez.api_key = self.apikey

    def _fetch_abstract_ids(self, disease: str, max_articles: int) -> List[str]:
        """
        Fetches PubMed IDs for articles related to the given disease.

        Args:
            disease (str): Disease keyword to search in PubMed.
            max_articles (int): Maximum number of article IDs to retrieve.

        Returns:
            List[str]: A list of PubMed article IDs.
        """
        try:
            logging.info(f"Searching PubMed for '{disease}'...")
            with Entrez.esearch(db="pubmed", term=disease, retmax=max_articles) as handle:
                record = Entrez.read(handle)
            ids = record.get("IdList", [])
            logging.info(f"Found {len(ids)} articles.")
            return ids
        except Exception as e:
            raise CustomException(error_message=e, error_details=sys)

    def _fetch_abstract_by_id(self, pmid: str) -> str:
        """
        Fetches the abstract text for a given PubMed ID.

        Args:
            pmid (str): PubMed ID of the article.

        Returns:
            str: Abstract text.
        """
        try:
            with Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text") as handle:
                return handle.read()
        except Exception as e:
            logging.warning(f"Failed to fetch abstract for PMID: {pmid}")
            return ""

    def _create_table(self, cursor) -> None:
        """
        Creates the 'Abstracts' table in the SQLite database if it doesn't exist.

        Args:
            cursor: SQLite database cursor.
        """
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS Abstracts (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                abstract TEXT NOT NULL
            )
        """)

    def _save_abstract(self, cursor, abstract: str) -> None:
        """
        Inserts a single cleaned abstract into the database.

        Args:
            cursor: SQLite database cursor.
            abstract (str): Cleaned abstract text.
        """
        cursor.execute("INSERT INTO Abstracts (abstract) VALUES (?)", (abstract,))

    def load_abstracts_to_db(self, disease: str, max_articles: int = 10) -> None:
        """
        Fetches abstracts for a disease from PubMed and stores them into a SQLite database.

        Args:
            disease (str): Disease keyword to search for.
            max_articles (int, optional): Number of articles to retrieve. Defaults to 10.
        """
        try:
            os.makedirs(self.data_folder_path, exist_ok=True)
            con = sqlite3.connect(self.db_path)
            cur = con.cursor()

            logging.info("Setting up database...")
            self._create_table(cur)

            ids = self._fetch_abstract_ids(disease, max_articles)
            logging.info("Fetching and storing abstracts...")

            for pmid in tqdm(ids):
                abstract = self._fetch_abstract_by_id(pmid)
                if abstract.strip():
                    cleaned_abstract = re.sub(r'\s+', ' ', abstract.replace('\n', ' ')).strip()
                    self._save_abstract(cur, cleaned_abstract)
                    logging.info(f"Stored abstract for PMID: {pmid}")
                else:
                    logging.warning(f"No abstract text for PMID: {pmid}")

            con.commit()
            logging.info("All abstracts saved successfully.")
        except Exception as e:
            raise CustomException(error_message=e, error_details=sys)
        finally:
            con.close()
            logging.info("Database connection closed.")

    def fetch_abstracts(self, disease: str, max_articles: int) -> List[str]:
        """
        Fetches and returns abstracts without storing them in the database.

        Args:
            disease (str): Disease keyword to search for.
            max_articles (int): Number of abstracts to fetch.

        Returns:
            List[str]: List of cleaned abstract texts.
        """
        try:
            ids = self._fetch_abstract_ids(disease, max_articles)
            abstracts = []

            logging.info("Fetching abstracts without DB storage...")
            for pmid in tqdm(ids):
                abstract = self._fetch_abstract_by_id(pmid)
                if abstract.strip():
                    cleaned_abstract = re.sub(r'\s+', ' ', abstract.replace('\n', ' ')).strip()
                    abstracts.append(cleaned_abstract)
            logging.info("Abstracts fetched successfully.")
            return abstracts
        except Exception as e:
            raise CustomException(error_message=e, error_details=sys)


if __name__ == "__main__":
    data_injection = DataInjection()
    data_injection.load_abstracts_to_db(disease="cancer", max_articles=100_000)
