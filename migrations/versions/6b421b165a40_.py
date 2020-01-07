"""empty message

Revision ID: 6b421b165a40
Revises: 35e0098d7400
Create Date: 2019-12-10 14:29:08.246466

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '6b421b165a40'
down_revision = '35e0098d7400'
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.add_column('gene_circuit', sa.Column('max_expression', sa.Integer(), nullable=True))
    op.add_column('gene_circuit', sa.Column('time_span', sa.Integer(), nullable=True))
    op.create_index(op.f('ix_gene_circuit_max_expression'), 'gene_circuit', ['max_expression'], unique=False)
    op.create_index(op.f('ix_gene_circuit_time_span'), 'gene_circuit', ['time_span'], unique=False)
    op.drop_index('ix_gene_circuit_circuit_name', table_name='gene_circuit')
    op.create_index(op.f('ix_gene_circuit_circuit_name'), 'gene_circuit', ['circuit_name'], unique=False)
    op.drop_column('gene_circuit', 'output_node')
    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.add_column('gene_circuit', sa.Column('output_node', sa.INTEGER(), nullable=True))
    op.drop_index(op.f('ix_gene_circuit_circuit_name'), table_name='gene_circuit')
    op.create_index('ix_gene_circuit_circuit_name', 'gene_circuit', ['circuit_name'], unique=1)
    op.drop_index(op.f('ix_gene_circuit_time_span'), table_name='gene_circuit')
    op.drop_index(op.f('ix_gene_circuit_max_expression'), table_name='gene_circuit')
    op.drop_column('gene_circuit', 'time_span')
    op.drop_column('gene_circuit', 'max_expression')
    # ### end Alembic commands ###
