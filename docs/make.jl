using SpeedyWeather
using WilliamsonTests
using Documenter

DocMeta.setdocmeta!(WilliamsonTests, :DocTestSetup, :(using WilliamsonTests); recursive=true)

makedocs(;
    modules=[WilliamsonTests],
    authors="Milan Klöwer <milankloewer@gmx.de> and contributors",
    sitename="WilliamsonTests.jl",
    format=Documenter.HTML(;
        canonical="https://speedyweather.github.io/WilliamsonTests.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Williamson's Tests" => [
            "test1.md",
        ],
    ],
)

deploydocs(
    repo="github.com/SpeedyWeather/WilliamsonTests.jl",
    devbranch="main",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
)
