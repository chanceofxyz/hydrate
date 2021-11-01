# Hydrate
Generate n simulation trials for a SIP in a SIPMath json file

CAUTION - This was hastily created for a hackathon. It worked with this test file
https://ipfs.io/ipfs/QmTn16U9YtbMeGkYWtRFZ5XrNA7tiKQYyRSgX4Es3TWouB

### INSTALL

TBD - just use the hydrate.js file for now. It has a dependency on jstat

```js
hydrateLibrary(LibraryIn, tokenIn, numberOfTrialsRequested)
```

...where LibraryIn is a json object for an entire SIPMath library, tokenIn is the name of a SIP in that library, and number of trials is well the numer of monte carlo trials you want returned.

### RUN

You can run `hydrateLibrary(LibraryIn, tokenIn, numberOfTrialsRequested)` by opening `run.html` in a web browser and examining console logs.
