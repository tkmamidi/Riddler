from Bio import Entrez
import os
from embedchain import App  # OpenSourceApp  # App -> for OpenAI API
from embedchain.config import ChatConfig


def get_pmc_ids(keywords):
    """
    Get a list of PMC IDs from PubMed Central (PMC) for a given search term.
    Parameters:
        keywords (str): Search term.
        n_articles (int): Number of articles to return.
    """
    Entrez.email = "your.email@example.com"
    # Search PMC for articles matching keywords
    handle = Entrez.esearch(db="pmc", term=keywords,
                            mindate=2005, maxdate=2023)
    record = Entrez.read(handle)
    print(f"Found {record['Count']} articles.")  # Print number of articles found

    # Get PMC IDs for articles matching keywords for a given number of articles
    n_articles = input("Enter number of articles to train from: ")
    handle = Entrez.esearch(db="pmc", sort="relevance", term=keywords, retmax=n_articles)
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
    chat_config = ChatConfig(stream=True)  # Set chat config to stream response

    # Get user input for keywords and number of articles to train from
    keywords = input("Enter your keywords: ")

    # Get PMC IDs from PubMed Central
    pmc_ids = get_pmc_ids(keywords)
    print(f"Selected {len(pmc_ids)} articles.")

    # Add articles to Riddle app to train from
    if pmc_ids:
        for pmc_id in pmc_ids:
            riddle.add("web_page", f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}")

    # Start Riddle app
    while True:
        question = input(f"\nRiddle me about (type 'end' to quit): ")
        if question == "end":
            break
        else:
            resp = riddle.chat(question, chat_config)
            answer = ""
            for chunk in resp:
                answer = answer + chunk
                print(chunk, end="", flush=True)
            print("\n")
            riddle.add_local("qna_pair", (question, answer))
