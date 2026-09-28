import asyncio, re, json, csv, time
from playwright.async_api import async_playwright
SP="/path/to/scratch"
RES="/path/to/revision/results/v4/"

async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch()
        ctx=await b.new_context(viewport={"width":1700,"height":1300})
        pg=await ctx.new_page()
        await pg.goto("https://compbio.cn/timer3/", wait_until="networkidle", timeout=300000)
        await pg.wait_for_timeout(8000)
        # Exploration -> Gene_Correlation
        await pg.click("a[data-value='Exploration']"); await pg.wait_for_timeout(2500)
        await pg.click("a[data-value='Gene_Correlation']"); await pg.wait_for_timeout(6000)
        await pg.screenshot(path=SP+"/timer_genecorr_0.png", full_page=False)
        vis=await pg.evaluate("""Array.from(document.querySelectorAll('#Gene_Correlation input, #Gene_Correlation select, #Gene_Correlation button, .selectize-input')).map(e=>{
            const r=e.getBoundingClientRect(); return {tag:e.tagName, id:e.id, cls:e.className, vis:(r.width>0&&r.height>0), html:e.outerHTML.slice(0,140)}}).filter(x=>x.vis)""")
        for v in vis[:30]: print("CTRL", v['tag'], v['id'], '|', v['html'][:110])
        await b.close()
asyncio.run(main())
