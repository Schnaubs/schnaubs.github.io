const express = require("express");
const axios = require("axios");
require("dotenv").config();

const app = express();
const PORT = 3000;

app.use(express.static("public")); // serve HTML

app.get("/api/card/:name", async (req, res) => {
  const cardName = req.params.name;

  try {
    const response = await axios.get(
      `https://api.pokemontcg.io/v2/cards?q=name:${cardName}`,
      {
        headers: {
          "X-Api-Key": process.env.POKEMON_API_KEY,
        },
      }
    );

    res.json(response.data.data);
  } catch (error) {
    res.status(500).json({ error: "Failed to fetch card data" });
  }
});

app.listen(PORT, () =>
  console.log(`Server running at http://localhost:${PORT}`)
);
