# Riddler

***!!! For research purposes only !!!***

Riddler is an LLM that looks up for free access research articles and answers any question you may have about a keyword.

## Usage

Riddler can be accessed at this streamlit [site](https://cgds-riddler.streamlit.app/).

## Installation

Installation simply requires fetching the source code. Following are required:

- Git

To fetch source code, change in to directory of your choice and run:

```sh
git clone https://github.com/uab-cgds-worthey/Riddler.git
```

### Set OpenAI API key

`export OPENAI_API_KEY="your_key_here"`

*Note*: No spaces between and use quotations. Get your key from [here](https://platform.openai.com/account/api-keys)


### Requirements

*OS:*

Currently works only in Mac OS. Docker versions may need to be explored later to make it useable in Mac (and
potentially Windows).

*Tools:*

- Python 3.9
- Pip3

*Environment:*

- [python virtual environment](https://docs.python.org/3/tutorial/venv.html)

### Install required packages

Change in to root directory and run the commands below:

```sh
# create environment. Needed only the first time.
pip3 install -r requirements.txt
```

## Steps to run

### Create an embedding model

Run the below command with 'keywords' as you'd search in pubmed and choose number of articles to train from.

`python src/riddler.py`

### Chat with Riddler

`streamlit run src/streamlit.py`

## Contact information

For issues, please send an email with clear description to


|Name | contact|
------|--------|
Tarun Mamidi | tmamidi@uab.edu
Shaurita Hutchins | sdhutchins@uab.edu

