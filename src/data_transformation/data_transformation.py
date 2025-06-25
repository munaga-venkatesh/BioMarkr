import os
import re
import sys
import nltk
import sqlite3 
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool
from nltk.corpus import stopwords
from nltk.stem import WordNetLemmatizer
from nltk.tokenize import word_tokenize
from src.utils.logger import logging
from src.utils.exception import CustomException


nltk.download('punkt')
nltk.download('stopwords')
nltk.download('wordnet')
nltk.download('punkt_tab')


class DataTransformation:
    def __init__(self, data_path: str) -> None:
        self.data_path = data_path
        self.entity_dict = {
            "breast cancer": "CANCER_TYPE",
            "lung cancer": "CANCER_TYPE",
            "colorectal cancer": "CANCER_TYPE",
            "tp53": "GENE",
            "brca1": "GENE",
            "her2": "BIOMARKER",
            "ca-125": "BIOMARKER",
            "chemotherapy": "TREATMENT",
            "immunotherapy": "TREATMENT",
            "tamoxifen": "DRUG",
            "cisplatin": "DRUG",
            "metastasis": "DISEASE",
            "braf v600e": "MUTATION",
            "kras g12d": "MUTATION"
        }

    def _extract_data_from_sql_to_df(self) -> pd.DataFrame:
        try:
            with sqlite3.connect(self.data_path) as con:
                if con:
                    try:
                        df = pd.read_sql_query("SELECT * FROM Abstracts", con, index_col="id")
                        return df
                    except pd.io.sql.DatabaseError:  # Table doesn't exist or is empty
                        logging.warning("Abstracts table not found or is empty.")
                        return pd.DataFrame()  # Return empty DataFrame
                else:
                    logging.warning("Database connection is None.")
                    return pd.DataFrame()
        except Exception as e:
            raise CustomException(error_details=sys, error_message=e)
    
    def _preprocess_single_token(self, abstract) -> list[str]:
        lemmatizer = WordNetLemmatizer()
        lem_text = lemmatizer.lemmatize(abstract)
        token_text = word_tokenize(lem_text)
        stop_words = set(stopwords.words('english'))  # Set for faster lookup
        token_text = [w for w in token_text if w.lower() not in stop_words]
        return token_text

    def _preprocess_abstracts_token(self) -> list:
        df = self._extract_data_from_sql_to_df()
        if df.empty:
            return []
        else:
            abstracts = df['abstract']
            with Pool() as pool:
                results = list(pool.map(self._preprocess_single_token, abstracts))
            return results
    
    
if __name__ == "__main__":
    data_path = os.path.join(os.getcwd(), "data", "abstracts.db")
    data_transfromation = DataTransformation(data_path= data_path)
    import time
    start_time = time.time()
    tokens = data_transfromation._preprocess_abstracts_token()
    end_time = time.time()
    print(f"⏱️ Time taken: {end_time - start_time:.2f} seconds")
    start_time = time.time()
    tokens = data_transfromation._preprocess_abstracts_token_for()
    end_time = time.time()
    print(f"⏱️ Time taken: {end_time - start_time:.2f} seconds")
