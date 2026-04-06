import os
import sys
import urllib.request
import re

try:
    from google import genai
    from google.genai import types
except ImportError:
    print("Please install the Google GenAI SDK: pip install google-genai")
    sys.exit(1)

# --- Tools for the AI ---

def read_local_file(filepath: str) -> str:
    """
    Reads the contents of a local file. Use this to inspect the user's code.
    
    Args:
        filepath: The path to the file to read.
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
            # Truncate if too large to save tokens
            if len(content) > 50000:
                return content[:50000] + "\n...[TRUNCATED TO SAVE TOKENS]..."
            return content
    except Exception as e:
        return f"Error reading file {filepath}: {e}"

def list_local_directory(path: str = ".") -> str:
    """
    Lists the files and directories in the given path.
    Use this to understand the project structure before reading files.
    
    Args:
        path: The directory path to list. Defaults to current directory.
    """
    try:
        items = os.listdir(path)
        return "\n".join(items)
    except Exception as e:
        return f"Error listing directory {path}: {e}"

def read_web_page(url: str) -> str:
    """
    Fetches the content of a web page as plain text. Use this to read PROS documentation.
    For PROS v4 docs, the base URL is https://pros.cs.purdue.edu/v5/pros-4/
    
    Args:
        url: The web address to fetch.
    """
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=10) as response:
            html = response.read().decode('utf-8')
            # Rudimentary HTML tag stripping to save tokens
            text = re.sub('<[^<]+?>', ' ', html)
            text = re.sub(r'\s+', ' ', text)
            if len(text) > 40000:
                return text[:40000] + "\n...[TRUNCATED]"
            return text
    except Exception as e:
        return f"Error fetching URL {url}: {e}"

def main():
    # 1. Check for API key
    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        print("Error: GEMINI_API_KEY environment variable not found.")
        print("\nPlease set it before running this script.")
        print("> PowerShell: $env:GEMINI_API_KEY=\"your_key_here\"")
        print("> CMD:        set GEMINI_API_KEY=your_key_here")
        sys.exit(1)

    print("Initializing PROS Assistant...")
    client = genai.Client(api_key=api_key)

    # 2. Define tools: Local File Access + Web Access (for PROS DOCS)
    tools = [
        read_local_file, 
        list_local_directory, 
        read_web_page
    ]

    # Gather project context (files and odom.hpp explicitly)
    print("Gathering project context to save API requests...", end="\r")
    project_files = []
    for root, _, files in os.walk("."):
        if any(ignored in root for ignored in [".git", ".vscode", "bin", "build"]):
            continue
        for file in files:
            if file.endswith(('.cpp', '.hpp', '.h', '.c')):
                project_files.append(os.path.normpath(os.path.join(root, file)))
    
    context_str = "Project C++ Files Overview:\n" + "\n".join(project_files) + "\n\n"
    
    # Pre-load odom.hpp if it exists (since this chatbot is for the particle filter)
    odom_loaded = False
    for pf in project_files:
        if "odom.hpp" in pf or "odom.h" in pf:
            context_str += f"--- PRE-LOADED CONTENTS OF {pf} ---\n"
            try:
                with open(pf, "r", encoding="utf-8") as f:
                    context_str += f.read()[:25000] # Truncate to keep prompt size reasonable
            except Exception as e:
                context_str += f"Could not read {pf}: {e}"
            context_str += "\n-----------------------------------------\n"
            odom_loaded = True
            break
            
    print("                                                 \r", end="")

    # 3. Define the instruction to keep token usage low and mandate PROS docs
    system_instruction = (
        "You are an expert C++ coding assistant specializing in the VEX PROS framework (v4/v5). "
        "Your primary goal is to help users set up and use the particle filter from their project.\n\n"
        f"Here is the local project context automatically loaded for you to save API calls:\n{context_str}\n"
        "1. You have tools `read_local_file` and `list_local_directory`. NEVER ask the user to paste code! "
        "If you need to see a file that isn't pre-loaded above, call `read_local_file` to read it yourself.\n"
        "2. To ensure accuracy with PROS APIs, use `read_web_page` for the official documentation. "
        "Use URLs starting with 'https://purduesigbots.github.io/pros-doxygen-docs/api.html'.\n"
        "3. Provide concise, accurate C++ code, and explain configurations clearly."
    )

    try:
        chat = client.chats.create(
            model="gemini-2.5-flash",
            config=types.GenerateContentConfig(
                tools=tools,
                system_instruction=system_instruction,
                temperature=0.2, # Low temperature for more factual coding responses
            )
        )
    except Exception as e:
        print(f"Error creating chat session: {e}")
        sys.exit(1)

    print("\n" + "="*60)
    print("🤖 PROS GenAI Chat Assistant Ready!")
    print(f"📂 Current Directory: {os.getcwd()}")
    print("Type 'quit' or 'exit' to stop.")
    print("="*60 + "\n")

    while True:
        try:
            user_input = input("You: ")
            if user_input.strip().lower() in ['quit', 'exit']:
                break
            if not user_input.strip():
                continue
                
            print("Gemini is thinking (investigating files/docs)...", end="\r")
            
            # Send message. The SDK natively handles function calls under the hood!
            response = chat.send_message(user_input)
            
            # Clear the loading message
            print("                                                  \r", end="")
            
            print("Gemini:\n" + response.text + "\n")
            
        except KeyboardInterrupt:
            print("\nExiting...")
            break
        except Exception as e:
            print(f"\nAn error occurred during generation: {e}")

if __name__ == "__main__":
    main()
