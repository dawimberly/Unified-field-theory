# Bringing your ChatGPT chats into the Unified-field-theory repo

This folder holds a small tool that turns your ChatGPT history into readable files.
Each chat becomes one Markdown page, for example `2025-03-02_hypersphere-and-e8.md`. An `INDEX.md` lists them all.
The pages are meant to go in the repo under `docs/chatgpt/`.

## Step 1: ask ChatGPT for your data
These steps come from OpenAI's help page "Exporting your ChatGPT history and data", checked on 24 Nov 2026.

1. Sign in at chatgpt.com, the same account you used for the physics chats.
2. Click your profile picture, then **Settings**.
3. Click **Data controls**.
4. Next to **Export data**, click **Export**, then **Confirm export**.
   - On the phone app the path is Settings → Data controls → Export data → Confirm export.
5. OpenAI emails (or texts) you when it's ready. That usually takes minutes, but it can take **up to 7 days**.
6. Open the email and click **Download data export** while signed in to the same account.
   - **The link stops working 24 hours after the email arrives**, so download it the same day.
   - If it expires, just request a new one.

You get a `.zip` file. Inside is `conversations.json` (or several `conversations-000.json`, `-001.json` files for big histories), plus any images you uploaded.
The zip can contain personal account details, so don't post it anywhere public.

Self-service export works for Free, Go, Plus and Pro accounts. If you're on a Business or Enterprise workspace, the workspace owner has to do it.

## Step 2: get the zip to me
Pick whichever is easier:

- **A. Attach it in our chat.** Drag the zip into the chat window and send it. If it's too big to attach, use option B.
- **B. Leave it on your PC (DeezDesktop) and tell me where.** Downloads usually land in `C:\Users\<your name>\Downloads\`. Send me the full path, e.g. `C:\Users\Dan\Downloads\chatgpt-export.zip`. I'll copy it over myself. You don't need to unzip it.

## Step 3: what I'll do with it
- Run the converter. By default it keeps only physics-related chats (keywords like physics, field, gravity, quantum, dimension, hypersphere, E8, unified, tesseract, brane…). If you'd rather keep everything, or add your own words, just tell me.
- Put the pages in `docs/chatgpt/` on a branch and show you the list before anything is pushed.
- Nothing goes into the repo until you say so.

## For reference: running it yourself
Python 3.8 or newer, no extra installs:

```
python convert_chatgpt_export.py chatgpt-export.zip -o docs/chatgpt --filter
python convert_chatgpt_export.py chatgpt-export.zip -o docs/chatgpt                 # keep every chat
python convert_chatgpt_export.py chatgpt-export.zip -o docs/chatgpt --filter E8 tesseract gravity
```

Extra switches:
- `--all-branches` also includes answers you regenerated or questions you edited.
- `--include-tools` shows code-run output.
- `--copy-assets` copies your uploaded images next to the pages.

Times in the files are in UTC. `python test_convert.py` checks the tool against a small made-up export.
