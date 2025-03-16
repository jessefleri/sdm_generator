# Species distribution models generator
 
R shiny application to create species distribution models from user data. Future iterations of this application will *hopefully* integrate a prompt style chatbot that will automate data scraping from existing APIs that serve ecological occurence data.

Link to project: https://jessefleri.shinyapps.io/sdm_generator/

How It's Made: `R` 

Core R packages: `shiny`, `ranger`, `leaflet`, `sf`, `geodata`

This is somewhere between a passion project and an opportunity to showcase what can be done in R. I'd love to see more folks consider R as a language capable of, dare I say, full stack engineering. 

Lessons learned:
- grepl pattern matching in to autopopulate columns selection
- Progress bars modals in the UI
- more conventional dashboard style than Ive built in the past 
