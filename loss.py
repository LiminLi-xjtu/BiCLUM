
import torch
import torch.nn as nn

def mask_correlated_samples(batch_size):
    N = 2 * batch_size
    mask = torch.ones((N, N))
    mask = mask.fill_diagonal_(0)
    for i in range(batch_size):
        mask[i, batch_size + i] = 0
        mask[batch_size + i, i] = 0
    mask = mask.bool()
    return mask


def contrastive_loss(batch_size, temperature, z_i, z_j):
    N = 2 * batch_size
    z = torch.cat((z_i, z_j), dim=0)

    mask = mask_correlated_samples(batch_size)

    sim = torch.matmul(z, z.T) / temperature
    sim_i_j = torch.diag(sim, batch_size)
    sim_j_i = torch.diag(sim, -batch_size)

    positive_samples = torch.cat((sim_i_j, sim_j_i), dim=0).reshape(N, 1)
    negative_samples = sim[mask].reshape(N, -1)

    labels = torch.zeros(N).to(positive_samples.device).long()
    logits = torch.cat((positive_samples, negative_samples), dim=1)
    criterion = nn.CrossEntropyLoss(reduction="sum")
    loss = criterion(logits, labels)
    loss /= N

    return loss


