const fs = require('fs')

/*
    * Reads all directory names from env variable SAMPLES_PATH
    * Returns an array of strings
    * Each string is a directory name
    */
async function getAllSamples() {
    const samplesPath = process.env.SAMPLES_PATH
    if (samplesPath === undefined) {
        throw new Error("Env variable SAMPLES_PATH is undefined")
    }

    const samples = await fs.promises.readdir(samplesPath)
    return samples
}

module.exports = getAllSamples

