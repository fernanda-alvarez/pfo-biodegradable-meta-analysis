import fs from "node:fs";
import path from "node:path";
import vm from "node:vm";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..");

const dataDir = path.join(root, "data");
const figDir = path.join(root, "outputs", "figures");

fs.mkdirSync(figDir, { recursive: true });

// -------------------------------------------------------------
// 1. Simple Offline Robvis
// -------------------------------------------------------------
const robvisHtml = fs.readFileSync(path.join(here, "Simple Offline Robvis", "index.html"), "utf8");
const rStart = robvisHtml.indexOf("<script>") + "<script>".length;
const rCut = robvisHtml.indexOf('\ndocument.getElementById("tool").addEventListener');
const rCore = robvisHtml.slice(rStart, rCut);

function createRobvisContext() {
  const ctx = { console };
  vm.createContext(ctx);
  vm.runInContext(rCore, ctx, { filename: "robvis-core.js" });
  return ctx;
}

// 1A. ROBINS-I
const robvisRobinsCtx = createRobvisContext();
robvisRobinsCtx.STATE.tool = "ROBINS-I v2";
const robinsCsv = fs.readFileSync(path.join(dataDir, "robvis_ROBINS-I_v2_permissive_tier2.csv"), "utf8");
const robinsImp = robvisRobinsCtx.importData(robinsCsv, "ROBINS-I v2");
console.log(`ROBINS-I parsed ${robinsImp.rows.length} rows`);

const robinsTfSvg = robvisRobinsCtx.trafficSVG(robinsImp.rows, "ROBINS-I v2");
const robinsBarSvg = robvisRobinsCtx.barSVG(robinsImp.rows, "ROBINS-I v2");

fs.writeFileSync(path.join(figDir, "robins_traffic_plot.svg"), robinsTfSvg);
fs.writeFileSync(path.join(figDir, "robins_summary_bar.svg"), robinsBarSvg);

// 1B. RoB 2
const robvisRob2Ctx = createRobvisContext();
robvisRob2Ctx.STATE.tool = "ROB2";
const rob2Csv = fs.readFileSync(path.join(dataDir, "robvis_RoB2_Zhang2026.csv"), "utf8");
const rob2Imp = robvisRob2Ctx.importData(rob2Csv, "ROB2");
console.log(`RoB 2 parsed ${rob2Imp.rows.length} rows`);

const rob2TfSvg = robvisRob2Ctx.trafficSVG(rob2Imp.rows, "ROB2");
const rob2BarSvg = robvisRob2Ctx.barSVG(rob2Imp.rows, "ROB2");

fs.writeFileSync(path.join(figDir, "rob2_traffic_plot.svg"), rob2TfSvg);
fs.writeFileSync(path.join(figDir, "rob2_summary_bar.svg"), rob2BarSvg);

// -------------------------------------------------------------
// 2. Simple JBI Appraisal Visualization
// -------------------------------------------------------------
const jbiHtml = fs.readFileSync(path.join(here, "Simple JBI Appraisal Visualization", "index.html"), "utf8");
const jStart = jbiHtml.indexOf("<script>") + "<script>".length;
const jCut = jbiHtml.indexOf("\n/* ================================================================ *\n *  DOM WIRING AND INTERACTIONS");
const jCore = jbiHtml.slice(jStart, jCut);

const jbiCtx = { console };
vm.createContext(jbiCtx);
vm.runInContext(jCore, jbiCtx, { filename: "jbi-core.js" });

const jbiCsv = fs.readFileSync(path.join(dataDir, "JBI_Case_Series_7_final_single_arm.csv"), "utf8");
const jbiImp = jbiCtx.parseCsvData(jbiCsv);
console.log(`JBI parsed ${jbiImp.rows.length} rows`);

const jbiTfSvg = jbiCtx.trafficSVG(jbiImp.rows);
const jbiBarSvg = jbiCtx.barSVG(jbiImp.rows);

fs.writeFileSync(path.join(figDir, "jbi_traffic_plot.svg"), jbiTfSvg);
fs.writeFileSync(path.join(figDir, "jbi_summary_bar.svg"), jbiBarSvg);

console.log("All risk of bias SVG figures generated successfully in outputs/figures/!");
