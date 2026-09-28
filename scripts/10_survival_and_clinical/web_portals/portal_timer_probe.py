import asyncio
from playwright.async_api import async_playwright
SP="/path/to/scratch"
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); ctx=await b.new_context(viewport={"width":1600,"height":1200})
        pg=await ctx.new_page()
        await pg.goto("https://compbio.cn/timer3/", wait_until="networkidle", timeout=240000)
        await pg.wait_for_timeout(6000)
        title=await pg.title(); print("TITLE:",title)
        tabs=await pg.evaluate("Array.from(document.querySelectorAll('a[data-value], .sidebar-menu a')).map(e=>(e.getAttribute('data-value')||'')+' :: '+e.textContent.trim()).slice(0,40)")
        for t in tabs: print("TAB:",t)
        await pg.screenshot(path=SP+"/timer_home.png", full_page=False)
        print("screenshot ok")
        await b.close()
asyncio.run(main())
