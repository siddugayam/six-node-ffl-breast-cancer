import asyncio, sys
from playwright.async_api import async_playwright
SP="/path/to/scratch"
GENE = sys.argv[1] if len(sys.argv) > 1 else "COL1A1"

async def sel(pg, sid, value):
    ctrl = "#%s + .selectize-control, #%s ~ .selectize-control" % (sid, sid)
    await pg.locator(ctrl).first.click(timeout=20000)
    await pg.wait_for_timeout(500)
    box = pg.locator(ctrl + " input").first
    await box.fill("")
    await box.type(value, delay=80)
    await pg.wait_for_timeout(2200)
    await pg.keyboard.press("Enter")
    await pg.wait_for_timeout(900)

async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch()
        pg = await (await b.new_context(viewport={"width": 1800, "height": 1500})).new_page()
        await pg.goto("https://compbio.cn/timer3/", wait_until="networkidle", timeout=300000)
        await pg.wait_for_timeout(8000)
        await pg.click("a[data-value='Immune']"); await pg.wait_for_timeout(2500)
        await pg.click("a[data-value='Immune_Gene']"); await pg.wait_for_timeout(9000)
        try: await sel(pg, "geneInput_gene", GENE)
        except Exception as e: print("gene sel:", repr(e)[:120])
        try: await sel(pg, "geneInput_infiltrates", "Cancer associated fibroblast")
        except Exception as e: print("infil sel:", repr(e)[:120])
        # purity adjustment on (default), cancer type ALL (default)
        try:
            await pg.get_by_role("button", name="Submit").first.click(timeout=20000)
            print("submitted")
        except Exception as e:
            print("submit failed:", repr(e)[:150])
        await pg.wait_for_timeout(70000)
        await pg.screenshot(path=SP + "/timer_immune_%s.png" % GENE, full_page=True)
        print("screenshot saved")
        txt = await pg.evaluate("document.body.innerText")
        i = txt.find("Instruction")
        print("TEXT AROUND RESULT:\n", txt[max(0, i-2500):i][:2500])
        await b.close()
asyncio.run(main())
