using PubChemReactions, HTTP, Test

# Local stand-in for PubChem: answers 503 while "throttled", like PubChem does under load.
function with_server(f, handler)
    server = HTTP.serve!(handler, "127.0.0.1", 0)
    return try
        f("http://127.0.0.1:$(HTTP.port(server))")
    finally
        close(server)
    end
end

@testset "pubchem_get outlasts a multi-second 503 throttling window" begin
    first_hit = Ref(NaN)
    handler = function (_)
        isnan(first_hit[]) && (first_hit[] = time())
        return time() - first_hit[] < 3.0 ? HTTP.Response(503, "busy") : HTTP.Response(200, "ok")
    end
    with_server(handler) do url
        res = PubChemReactions.pubchem_get(url)
        @test res.status == 200
        @test String(res.body) == "ok"
    end
end

@testset "pubchem_get does not retry a definitive 404" begin
    hits = Ref(0)
    with_server(_ -> (hits[] += 1; HTTP.Response(404, "not found"))) do url
        @test_throws HTTP.StatusError PubChemReactions.pubchem_get(url)
        @test hits[] == 1
    end
end
