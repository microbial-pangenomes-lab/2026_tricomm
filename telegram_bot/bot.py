import os
import asyncio
import logging

import httpx
from telegram import Update
from telegram.ext import Application, CommandHandler, ContextTypes

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)

TARGET_URL = "http://192.168.7.2:5000/"
CHECK_INTERVAL = 1800  # seconds

# chat_id -> asyncio.Task
monitor_tasks: dict[int, asyncio.Task] = {}


async def check_url() -> bool:
    try:
        async with httpx.AsyncClient() as client:
            resp = await client.get(TARGET_URL, timeout=10)
            return resp.status_code < 500
    except (httpx.RequestError, httpx.TimeoutException):
        return False


async def start(update: Update, context: ContextTypes.DEFAULT_TYPE) -> None:
    await update.message.reply_text(
        "Hi! I can monitor whether the chibio interface is reachable.\n\n"
        "Commands:\n"
        "/check — check the connection to the interface right now\n"
        f"/monitor — check every {CHECK_INTERVAL // 60} min, alert when down\n"
        "/stop — stop monitoring"
    )


async def check(update: Update, context: ContextTypes.DEFAULT_TYPE) -> None:
    up = await check_url()
    await update.message.reply_text(
        "✅ The interface is up!" if up else "❌ The interface is down!"
    )


async def _monitor_loop(chat_id: int, app: Application) -> None:
    while True:
        up = await check_url()
        if not up:
            await app.bot.send_message(chat_id, "⚠ The interface is down!")
        await asyncio.sleep(CHECK_INTERVAL)


async def monitor(update: Update, context: ContextTypes.DEFAULT_TYPE) -> None:
    chat_id = update.effective_chat.id
    if chat_id in monitor_tasks:
        await update.message.reply_text("Monitoring is already running.")
        return
    task = asyncio.create_task(_monitor_loop(chat_id, context.application))
    monitor_tasks[chat_id] = task
    await update.message.reply_text(
        f"Monitoring started. I'll check every {CHECK_INTERVAL // 60} minutes "
        "and alert you if the URL is down."
    )


async def stop(update: Update, context: ContextTypes.DEFAULT_TYPE) -> None:
    chat_id = update.effective_chat.id
    task = monitor_tasks.pop(chat_id, None)
    if task is None:
        await update.message.reply_text("Monitoring is not running.")
        return
    task.cancel()
    await update.message.reply_text("Monitoring stopped.")


def main() -> None:
    token = os.environ["TELEGRAM_BOT_TOKEN"]
    app = Application.builder().token(token).build()
    app.add_handler(CommandHandler("start", start))
    app.add_handler(CommandHandler("check", check))
    app.add_handler(CommandHandler("monitor", monitor))
    app.add_handler(CommandHandler("stop", stop))
    app.run_polling()


if __name__ == "__main__":
    main()
