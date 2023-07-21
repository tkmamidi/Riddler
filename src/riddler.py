from Bio import Entrez
import os
from embedchain import App


def get_pmc_ids(keywords, n_articles=10):
    """
    Get a list of PMC IDs from PubMed Central (PMC) for a given search term.
    Parameters:
        keywords (str): Search term.
        n_articles (int): Number of articles to return.
    """
    Entrez.email = "your.email@example.com"
    handle = Entrez.esearch(db="pmc", term=keywords, retmax=n_articles)
    record = Entrez.read(handle)
    id_list = record["IdList"]
    return id_list


def get_api_key():
    api_key = os.environ.get("OPENAI_API_KEY")
    if api_key is None:
        raise ValueError("OPENAI_API_KEY environment variable not set.")
    return None


if __name__ == "__main__":
    get_api_key()  # Check if API key is set
    riddle = App()  # Initialize Riddle app

    # Get user input for keywords and number of articles to train from
    keywords = input("Enter your keywords: ")
    n_articles = input("Enter number of articles to train from: ")

    # Get PMC IDs from PubMed Central
    pmc_ids = get_pmc_ids(keywords, n_articles)

    # Add articles to Riddle app to train from
    if pmc_ids:
        for pmc_id in pmc_ids:
            riddle.add("web_page", f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}")

    # Start Riddle app
    while True:
        question = input(f"\nAsk me something about {keywords} (press cmd+c to quit): ")
        print("\n" + riddle.query(question))
