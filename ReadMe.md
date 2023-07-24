# Riddler

***!!! For research purposes only !!!***

Riddler is a ChatGPT-like app that incorporates the capabilities of Large Language Models (LLMs) in combination with [embedchain](https://github.com/embedchain/embedchain) to answer any questions based on PubMed keywords.

Riddler has been trained on 1000 open-access PubMed Central articles using "cystic fibrosis" as the keyword.

## Usage

Riddler can be accessed at this streamlit [site](https://cgds-riddler.streamlit.app/).

## Installation

Installation simply requires fetching the source code. The following software is required:

- [git](https://git-scm.com/downloads)

To fetch the source code, change into the directory of your choice and run:

```sh
git clone https://github.com/uab-cgds-worthey/Riddler.git
cd Riddler/
```

### Set OpenAI API key

`export OPENAI_API_KEY="your_key_here"`

*Note*: No spaces between and use quotations. Get your key from [here](https://platform.openai.com/account/api-keys).

### Requirements

*OS:*

Currently works only in Mac OS. Docker versions may need to be explored later to make it usable on Mac (and
potentially Windows).

*Tools:*

- Python 3.9
- pip3

*Environment:*

- [python virtual environment](https://docs.python.org/3/tutorial/venv.html)

### Create environment & install required packages

Change in to root directory and run the commands below:

```bash
# Create an environment. Needed only the first time.
python3 -m venv riddler-env
source riddler-env/bin/activate
pip3 install -r requirements.txt
```

## Steps to run

### Create an embedding model

Run the below command with 'keywords' as you'd search in pubmed and choose a number of articles to train from.

`python src/riddler.py`

### Chat with Riddler

`streamlit run src/streamlit.py`

## Authors

For issues, please send an email with a clear description to one of the authors.

|Name | contact|
------|--------|
Tarun Mamidi | <tmamidi@uab.edu>
Shaurita Hutchins | <sdhutchins@uab.edu>
