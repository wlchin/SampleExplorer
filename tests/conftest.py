import pytest

from sample_explorer.sample_explorer import Query_DB
from tests.common import SEMANTIC_DB, TRANSCRIPTOME_DB, build_mini_archs4_h5


@pytest.fixture(scope="session")
def archs4_h5(tmp_path_factory):
    return str(build_mini_archs4_h5(tmp_path_factory.mktemp("archs4") / "mini_archs4.h5"))


@pytest.fixture(scope="session")
def query_db():
    return Query_DB(SEMANTIC_DB, TRANSCRIPTOME_DB)


@pytest.fixture(scope="session")
def query_db_archs4(archs4_h5):
    return Query_DB(SEMANTIC_DB, TRANSCRIPTOME_DB, archs4_h5)
