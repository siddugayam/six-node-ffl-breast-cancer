import asyncio, re
from playwright.async_api import async_playwright
SP="/path/to/scratch"
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(args=["--ignore-certificate-errors"])
        ctx=await b.new_context(ignore_https_errors=True, viewport={"width":1500,"height":1100})
        pg=await ctx.new_page()
        bodies=[]
        pg.on("request", lambda r: bodies.append((r.url,r.post_data)) if r.method=="POST" else None)
        await pg.goto("https://bcgenex.ico.unicancer.fr/BC-GEM/GEM-Requete.php?mode=1", wait_until="networkidle", timeout=180000)
        oc=await pg.evaluate("document.getElementById('Verifier').getAttribute('onclick')")
        print("VERIFIER ONCLICK:", oc)
        await pg.click("#token-input-Gene2TextID2")
        await pg.keyboard.type("PTEN", delay=120)
        await pg.wait_for_timeout(3500)
        drop=await pg.evaluate("Array.from(document.querySelectorAll('.token-input-dropdown li, ul li')).map(e=>e.textContent.trim()).filter(t=>t&&t.length<40).slice(0,15)")
        print("dropdown:",drop)
        try:
            await pg.click(".token-input-dropdown li:first-child", timeout=8000); print("clicked dropdown")
        except Exception as e:
            print("dropdown click failed:",e); await pg.keyboard.press("Enter")
        await pg.wait_for_timeout(1200)
        print("tokens:", await pg.evaluate("Array.from(document.querySelectorAll('.token-input-token')).map(e=>e.textContent.trim())"))
        await pg.click("#Verifier")
        await pg.wait_for_load_state("networkidle", timeout=300000)
        await pg.wait_for_timeout(4000)
        out=await pg.content(); open(SP+"/bcgem_result.html","w").write(out)
        print("URL:",pg.url,"len",len(out))
        await pg.screenshot(path=SP+"/bcgem_result.png", full_page=True)
        for u,d in bodies:
            print("POST",u); 
            if d: open(SP+"/bcgem_postbody.txt","w").write(d); print(d[:1500])
        await b.close()
asyncio.run(main())
