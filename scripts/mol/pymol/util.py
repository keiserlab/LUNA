
def mybio_entity_2selection(entity):
    params = {}
    if entity.level == 'S' or entity.level == 'M':
        params[''] = 'all'
    else:
        if entity.level == 'C':
            params['chain'] = entity.id
        else:
            if entity.level == 'R':
                params['resn'] = entity.resname
                params['res'] = str(entity.id[1]) + entity.id[2].strip()
                params['chain'] = entity.get_parent().id
            else:
                if entity.level == 'A':
                    residue = entity.get_parent()
                    params['id'] = entity.serial_number
                    params['name'] = entity.name
                    params['resn'] = residue.resname
                    params['res'] = str(residue.id[1]) + residue.id[2].strip()
                    params['chain'] = residue.get_parent().id
            return (' AND ').join(['%s %s' % (k, str(v)) for k, v in params.items()])
