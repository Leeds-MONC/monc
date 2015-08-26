package com.epcc.monc;

import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

public class HtmlProcessor {
	private static final String HTML_PREAMBLE="<html style=\"font-family: helvetica !important;\"><head><title>MONC unit test results</title>"+
			"<link href=\"jquery-ui.css\" rel=\"stylesheet\"><script src=\"jquery-1.10.2.js\"></script>"+
			"<script src=\"jquery-ui-1.10.4.custom.js\"></script><script>$(function() {$( \"#dialog\" ).dialog({autoOpen: false,width:1000, height:600});});"+
			"var displayInfoDialog = function(message, title) {  $(\"div#dialog\").html(message);$(\"#dialog\").dialog('option', 'title', title);$(\"#dialog\").dialog(\"open\");};</script>"+
			"</head><body><div id='dialog' ></div><h1>MONC Unit test results</h1>";
			
	private static final String HTML_POSTAMBLE="</body></html>";
	public static String generate(String timestamp, Map<String, CompilerTests> compilerTests, Set<String> unitTestCategories, 
			Map<String, String> unitTestURLs, String version, boolean overallPass, int totalCases, int totalPasses) {		
		StringBuilder sb = new StringBuilder();
		sb.append(HTML_PREAMBLE);
		sb.append("<h2>Revision number <a href='https://puma.nerc.ac.uk/trac/MONC/timeline' target='blank'style=\"text-decoration: none\">").append(version).append(
				"</a> <i><a href='http://www2.epcc.ed.ac.uk/~nbrown23/monc/docs/html/index.html' target='blank'style=\"text-decoration: none\">(docs)</a></i>,  run at ").append(timestamp).append("</h2>");
		sb.append("<div><b>Overall test run: <i>").append(String.valueOf(totalPasses)).append("/").append(String.valueOf(totalCases)).append(
				"</i></b><img width='15px' height='15px' src='");
		if (overallPass) {
			sb.append("http://upload.wikimedia.org/wikipedia/en/e/e4/Green_tick.png");
		} else {
			sb.append("http://img2.wikia.nocookie.net/__cb20130801192632/papersplease/images/f/fe/Red_cross_tick.png");
		}
		sb.append("'>");
		for (Entry<String, CompilerTests> entry : compilerTests.entrySet()) {
			sb.append("&nbsp;&nbsp;&nbsp;&nbsp;").append(entry.getKey()).append(": <i>").append(entry.getValue().getTotalTestCasePasses()).append("/").append(entry.getValue().getTotalTestCases()).append(
					"</i><img width='15px' height='15px' src='");
			if (entry.getValue().getPass()) {
				sb.append("http://upload.wikimedia.org/wikipedia/en/e/e4/Green_tick.png");
			} else {
				sb.append("http://img2.wikia.nocookie.net/__cb20130801192632/papersplease/images/f/fe/Red_cross_tick.png");
			}
			sb.append("'>");
		}
		sb.append("</div><br><table border=1><tr><td></td>");
		for (String specificCompiler : compilerTests.keySet()) {
			sb.append("<td style=\"padding-left:5px;padding-right:5px;\">").append(specificCompiler).append("</td>");
		}
		sb.append("</tr>");
		
		sb.append("<tr><td style=\"padding-left:5px;padding-right:5px;\">MONC build</td>");
		for (CompilerTests cT : compilerTests.values()) {
			boolean passed = cT.getPassCompile();
			String messages = cT.getCompileMessage();
			buildPassFailItem(sb, "MONC build", cT, passed, messages.replace("'", "`").replace("\"", "`"));
		}
		sb.append("</tr>");
		for (String specificTest : unitTestCategories) {
			sb.append("<tr><td style=\"padding-left:5px;padding-right:5px;\">");
			if (unitTestURLs.containsKey(specificTest)) sb.append("<a href='https://").append(unitTestURLs.get(specificTest)).append("' target='blank'>");
			sb.append(specificTest);
			if (unitTestURLs.containsKey(specificTest)) sb.append("</a>");
			sb.append("</td>");
			for (CompilerTests cT : compilerTests.values()) {
				boolean passed = cT.getTestSuccess(specificTest);
				String messages = cT.getTestMessage(specificTest);
				buildPassFailItem(sb, specificTest, cT, passed, messages);
			}
			sb.append("</tr>");
		}		
		sb.append("</table>");
		sb.append(HTML_POSTAMBLE);
		return sb.toString();
	}
	private static void buildPassFailItem(StringBuilder sb,
			String specificTest, CompilerTests cT, boolean passed,
			String messages) {
		sb.append("<td style=\"padding-left:5px;padding-right:5px;\" bgcolor='#").append(passed ? "00ff00" : "ff0000").append("'>");
		if (messages != null) sb.append("<a href=\"javascript:displayInfoDialog('"+messages+"', '"+specificTest+" ("+cT.getName()+") unit test details');\">");
		if (cT.getTestCases(specificTest) > 0) {
			sb.append(cT.getTestCasePasses(specificTest)).append("/").append(cT.getTestCases(specificTest));
		} else {
			sb.append(passed ? "Pass" : "Fail");
		}
		if (messages != null) sb.append("</a>");
		sb.append("</td>");
	}
}
