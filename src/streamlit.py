import streamlit as st
from embedchain import App  # OpenSourceApp  # App -> for OpenAI API
from embedchain.config import ChatConfig
import time

# Config the whole app
st.set_page_config(
    page_title="RIDDLER",
    page_icon="🧊",
    layout="wide",  # initial_sidebar_state="expanded",
)
st.title("Riddler")


@st.cache_resource
def initialize_app():
    riddle = App()  # Initialize Riddle app
    chat_config = ChatConfig(stream=True)  # Set chat config to stream response
    return riddle, chat_config


@st.cache_data
def sidebar_links():
    software_link_dict = {
        "GitHub Page": "https://github.com/uab-cgds-worthey/Riddler",
        "Entrez esearch": "https://biopython.org/docs/1.75/api/Bio.Entrez.html#Bio.Entrez.esearch",
        "embedchain": "https://github.com/embedchain/embedchain",
        "Streamlit": "https://streamlit.io",
    }

    st.sidebar.markdown("## Related Links")

    for link_text, link_url in software_link_dict.items():
        st.sidebar.markdown(f"[{link_text}]({link_url})")

    st.sidebar.markdown("---")

    st.sidebar.markdown(
        f"""
    ### Authors
    Please feel free to contact us with any issues, comments, or questions.""",
    )

    st.sidebar.markdown(
        f"""
    ##### Tarun Mamidi [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social&label=Follow%20%40TarunMamidi7)](https://twitter.com/TarunMamidi7)
    - Email:  <tmamidi@uab.edu>
    - GitHub: https://github.com/tkmamidi
    ##### Shaurita Hutchins [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social&label=Follow%20%40shauritacodes)](https://twitter.com/shauritacodes)
    - Email: <sdhutchins@uab.edu>
    - GitHub: https://github.com/sdhutchins
    ##### Dr. Liz Worthey [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social&label=Follow%20%40lizworthey)](https://twitter.com/lizworthey)
    - Email: <lworthey@uab.edu>
    - GitHub: https://github.com/uab-cgds-worthey
    """,
    )

    st.sidebar.markdown("---")

    st.sidebar.markdown("Copyright (c) 2023 CGDS")
    return None


riddle, chat_config = initialize_app()
# Initialize chat history
if "messages" not in st.session_state:
    st.session_state.messages = []

# Display chat messages from history on app rerun
for message in st.session_state.messages:
    with st.chat_message(message["role"]):
        st.markdown(message["content"])

# React to user input
if user_input := st.chat_input("Riddle me about ..."):
    # Add user message to chat history
    st.session_state.messages.append({"role": "user", "content": user_input})
    # Display user message in chat message container
    with st.chat_message("user"):
        st.markdown(user_input)
    # Display assistant response in chat message container
    with st.chat_message("assistant"):
        message_placeholder = st.empty()
        full_response = ""
    with st.chat_message("assistant"):
        message_placeholder = st.empty()
        full_response = ""

        with st.spinner("Riddler unriddling riddle ..."):
            assistant_response = riddle.chat(user_input, chat_config)
        # Simulate stream of response with milliseconds delay
        for chunk in assistant_response:
            full_response += chunk + ""
            time.sleep(0.05)
            # Add a blinking cursor to simulate typing
            message_placeholder.markdown(full_response + "▌")
        message_placeholder.markdown(full_response)
    st.session_state.messages.append({"role": "assistant", "content": full_response})


sidebar_links()
